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
      attr_reader :mutex
      attr_reader :read_queue, :data_queue
      attr_reader :process_data

      r_attribute :n_cpu, :int, "Number of CPUs"
      r_attribute :min_io, :int, "Minimum I/O tasks"
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
        @min_io = 2
        @n_cpu = 4
        @workers = []
        @mutex = Mutex.new
        @read_queue = ConcurrentLinkedQueue.new
        @data_queue = ArrayDeque.new
        @blocked_stack = ArrayDeque.new
        @hi_wat = 16
        @lo_wat = 4
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
          && (! in_overflow? && ! blocked_stack.is_empty) \
          && ! read_queue.is_empty
          blocked_stack.remove.unblock
        end
      end

      def in_overflow?
        if @overflow
          # already in overflow, have we reached low watermark?
          @overflow = @data_queue.size <= @lo_wat
        else
          # not in overflow, have we reached high watermark?
          @overflow = @data_queue.size >= @hi_wat
        end
        @overflow
      end

      def needs_readers?
        (n_reading < min_io \
         || n_reading < n_processing) \
        && ! read_queue.is_empty
      end

    end

    class Worker
      attr_reader :controller, :thread, :fd
      attr_accessor :should_run

      def initialize(controller)
        @controller = controller
        @fd = File.open(controller.path)
        @should_run = true
        @state = :read
        @new_state = @state
        @block_barrier = CyclicBarrier.new(2)
      end

      def start
        @thread = Thread.new { run }
      end

      def run
        begin
          controller.mutex.synchronize do
            controller.count_by_state[@state] += 1
          end
          begin
            $stderr.puts "Worker running, state = #{@state}."
            while should_run
              case @state
              when :read
                do_read
              when :block
                do_block
              when :process
                do_process
              else
                raise "invalid state: #{@state}"
              end
              if @new_state != @state
                controller.mutex.synchronize do
                  controller.count_by_state[@state] -= 1
                  controller.count_by_state[@new_state] += 1
                end
                # $stderr.puts "Worker: #{@state} -> #{@new_state}"
                @state = @new_state
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

      def unblock
        @new_state = :read
        @block_barrier.await
      end

      def do_block
        @block_barrier.reset
        controller.mutex.synchronize do
          controller.blocked_stack.add_first(self)
        end
        @block_barrier.await
      end

      def read_data(req)
        fd.seek(req.offset)
        req.data = fd.read(req.length)
      end

      def do_read
        req = controller.read_queue.poll
        read_data(req)
        controller.mutex.synchronize do
          if controller.n_processing < controller.n_cpu
            @new_state = :process
            @cur_data = req
            controller.unblock_if_appropriate
            # keep block
          elsif controller.needs_readers?
            # stay in :read
            controller.data_queue.add(req)
          else
            @new_state = :block
            controller.data_queue.add(req)
          end
        end
      end

      def get_data
        if @cur_data
          data = @cur_data
          @cur_data = nil
        else
          controller.mutex.synchronize do
            data = controller.data_queue.poll
          end
        end
        data
      end

      def do_process
        data = get_data
        controller.process_data.call(data)
        controller.mutex.synchronize do
          if (! controller.data_queue.isEmpty()) \
            && (controller.n_processing < controller.n_cpu)
            # stay in :process
            controller.unblock_if_appropriate
          elsif controller.needs_readers?
            @new_state = :read
          else
            @new_state = :block
          end
        end
      end

    end

  end
end
