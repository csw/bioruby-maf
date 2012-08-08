require 'set'

module Bio::MAF
  
  module JobRunner
    def JobRunner.create(n_parallel)
      if RUBY_PLATFORM == 'java'
        JThreadRunner.new(n_parallel)
      else
        ForkRunner.new(n_parallel)
      end
    end
  end

  class ForkRunner

    def initialize(n_parallel)
      @n_parallel = n_parallel
      @jobs = []
      @kids = Set.new
    end

    def add(&proc)
      @jobs << proc
    end

    def run
      until @jobs.empty? && @kids.empty?
        while can_start?
          start_job
        end
        await
      end
    end

    private
    
    def can_start?
      (! @jobs.empty?) && @kids.size < @n_parallel
    end

    def start_job
      job = @jobs.shift
      pid = fork()
      if pid
        # parent
        @kids << pid
      else
        # child
        begin
          job.call()
          exit 0
        rescue SystemExit
          raise
        rescue Exception
          LOG.error $!
          exit 1
        end
      end
    end

    def await
      pid = Process.wait
      unless @kids.delete?(pid)
        raise "Completion of unexpected job #{pid}!"  
      end
      if ! $?.success?
        raise "Job #{pid} failed with status #{status.exitstatus}!"
      end
    end

  end

  class JThreadRunner

    def initialize(n_parallel)
      @n_parallel = n_parallel
      @exec = java.util.concurrent.Executors.newFixedThreadPool(n_parallel)
      @ecs = java.util.concurrent.ExecutorCompletionService.new(@exec)
      @n = 0
    end

    def add(&blk)
      @ecs.submit(&blk)
      @n += 1
    end

    def run
      seen = 0
      until seen == @n
        f = @ecs.take()
        begin
          f.get()
        rescue Exception => e
          raise "Job failed: #{e}"
        end
        seen += 1
      end
      @exec.shutdown()
    end

  end

end
