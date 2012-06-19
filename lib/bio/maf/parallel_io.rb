require 'java'
require 'jmx4r'

module Bio::MAF
  module ParallelIO
    java_import java.util.concurrent.ConcurrentLinkedQueue
    java_import java.util.concurrent.CyclicBarrier
    java_import java.util.ArrayDeque
    java_import java.lang.management.ManagementFactory
    java_import javax.management.ObjectName

    class Controller < JMX::DynamicMBean
      r_attribute :path, :string, "File path"
      attr_reader :mutex, :workers
      attr_reader :read_queue, :data_queue
      attr_reader :process_data

      rw_attribute :n_cpu, :int, "Number of CPUs"
      rw_attribute :min_io, :int, "Minimum I/O tasks"

      r_attribute :n_workers, :int, "Number of worker threads"
      r_attribute :n_processing, :int, "Current number of processing threads"
      r_attribute :n_reading, :int, "Current number of reading threads"
      r_attribute :n_blocked, :int, "Current number of blocked threads"

      attr_reader :count_by_state
      attr_reader :blocked_stack

      def initialize(path)
        super()
        @path = path
        @n_workers = 12
        @min_io = 1
        @n_cpu = 4
        @workers = []
        @mutex = Mutex.new
        @read_queue = ArrayDeque.new
        @data_queue = ArrayDeque.new
        @blocked_stack = ArrayDeque.new
        @hi_wat = 256
        @lo_wat = 64
        @overflow = false
        @count_by_state = {
          :read => 0,
          :block => 0,
          :process => 0
        }
        #oname = ObjectName.new("ParallelIO:path=#{path}")
        oname = ObjectName.new("ParallelIO:a=b")
        $stderr.puts "registering #{self.getMBeanInfo().getDescription()} as #{oname}"
        ManagementFactory.platform_mbean_server.register_mbean(self, oname)
      end

      def on_data(&block)
        @process_data = block
      end

      def start
        @n_workers.times do
          w = Worker.new(self)
          @workers << w
          w.start
        end
        $stderr.puts "started #{@n_workers} workers."
      end

      def finished?
        mutex.synchronize do
          read_queue.is_empty \
          && data_queue.is_empty \
          && count_by_state[:block] == n_workers
        end
      end

      def shutdown
        mutex.synchronize do
          until blocked_stack.is_empty
            w = blocked_stack.remove
            w.should_run = false
            w.unblock
          end
        end
      end

      def dump_state
        s = ''
        mutex.synchronize do
          [:read, :process, :block].each do |state|
            s << sprintf("%s: %2d  ", state.to_s, count_by_state[state])
          end
          s << sprintf("data queue: %3d ", data_queue.size)
          s << "[O] " if in_overflow?
          if read_queue.is_empty()
            s << "read queue empty"
          else
            s << "read queue not empty"
          end
        end
        $stderr.puts s
      end

      def n_processing
        count_by_state[:process]
      end
      def jmx_get_n_processing
        javax.management.Attribute.new("n_processing",
                                       n_processing)
      end

      def n_reading
        count_by_state[:read]
      end
      def jmx_get_n_reading
        javax.management.Attribute.new("n_reading",
                                       n_reading)
      end

      def n_blocked
        count_by_state[:block]
      end
      def jmx_get_n_blocked
        javax.management.Attribute.new("n_blocked",
                                       n_blocked)
      end

      def unblock_if_appropriate
        # caller must hold mutex
        if ((n_reading < min_io) || (n_reading < n_processing)) \
          && ! in_overflow? \
          && ! blocked_stack.is_empty \
          && ! read_queue.is_empty

          blocked_stack.remove.unblock
          # $stderr.puts "unblocking a thread"
        end
      end

      def in_overflow?
        if @overflow
          # already in overflow, have we reached low watermark?
          @overflow = (@data_queue.size >= @lo_wat)
        else
          # not in overflow, have we reached high watermark?
          @overflow = (@data_queue.size >= @hi_wat)
        end
        @overflow
      end

      def needs_readers?
        (n_reading < min_io || n_reading < n_processing) \
        && (! read_queue.is_empty)
        #&& (! in_overflow?) \
      end

    end

    class Worker
      attr_reader :controller, :thread, :fd
      attr_accessor :should_run

      def initialize(controller)
        @controller = controller
        @fd = File.open(controller.path)
        @should_run = true
        @state = :start
        @new_state = @state
        @block_barrier = CyclicBarrier.new(2)
      end

      def start
        @thread = Thread.new { run }
      end

      def transition_to(new_state)
        unless @cur_data
          case new_state
          when :read
            @cur_data = controller.read_queue.remove
          when :process
            @cur_data = controller.data_queue.remove
          end
        end

        if new_state != @state
          unless @state == :start
            controller.count_by_state[@state] -= 1
          end
          controller.count_by_state[new_state] += 1
          @state = new_state
        end
      end

      def run
        begin
          begin
            $stderr.puts "Worker running, state = #{@state}."
            while should_run
              case @state
              when :read
                do_read
              when :process
                do_process
              when :block
                do_block
              when :start
                do_start
              else
                raise "invalid state: #{@state}"
              end
            end
          ensure
            controller.mutex.synchronize do
              controller.count_by_state[@state] -= 1
            end
          end
        rescue
          $stderr.puts "Worker exiting: #{$!.class}: #{$!}"
          $stderr.puts $!.backtrace.join("\n")
        end
      end

      def do_start
        controller.mutex.synchronize do
          transition_to(:read)
          # if controller.needs_readers?
          #   transition_to(:read)
          # elsif ! controller.data_queue.isEmpty()
          #   transition_to(:process)
          # else
          #   transition_to(:block)
          # end
        end
      end

      def cur_req
        req = @cur_data
        unless req
          raise "No current request!"
        end
        @cur_data = nil
        req
      end

      def read_data(req)
        fd.seek(req.offset)
        req.data = fd.read(req.length)
      end

      def do_read
        req = cur_req()
        read_data(req)
        controller.mutex.synchronize do
          if controller.n_processing < controller.n_cpu
            # keep block and process it
            @cur_data = req
            transition_to(:process)
            controller.unblock_if_appropriate
          elsif controller.needs_readers?
            # stay in :read
            controller.data_queue.add(req)
            transition_to(:read)
            controller.unblock_if_appropriate
          else
            controller.data_queue.add(req)
            # $stderr.puts "read -> block"
            transition_to(:block)
          end
        end
      end

      def do_process
        data = cur_req()
        controller.process_data.call(data)
        controller.mutex.synchronize do
          if (! controller.data_queue.isEmpty()) \
            && (controller.n_processing < controller.n_cpu)
            transition_to(:process)
            controller.unblock_if_appropriate
          elsif controller.needs_readers?
            transition_to(:read)
            controller.unblock_if_appropriate
          else
            # $stderr.puts "process -> block"
            transition_to(:block)
          end
        end
      end

      def unblock
        transition_to(:read) if should_run
        @block_barrier.await
      end

      def do_block
        start = Time.now
        @block_barrier.reset
        controller.mutex.synchronize do
          controller.blocked_stack.add_first(self)
        end
        @block_barrier.await
        elapsed = Time.now - start
        # $stderr.printf("Unblocked after %.1fs.\n", elapsed)
      end

    end

  end
end
