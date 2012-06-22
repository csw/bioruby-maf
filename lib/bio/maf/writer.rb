module Bio::MAF

  class Writer
    attr_reader :f, :path

    def initialize(path)
      @path = path
      @f = File.open(path, 'w')
    end

    def flatten_vars(vars)
      vars.to_a.collect {|k, v| "#{k}=#{v}"}.join(" ")
    end

    def write_header(header)
      f.puts "##maf #{flatten_vars(header.vars)}"
      f.puts "##{header.alignment_params}" if header.alignment_params
    end

    def write_blocks(blocks)
      blocks.each do |block|
        write_block(block)
      end
      f.flush
    end

    def write_block(block)
      f.puts "a #{flatten_vars(block.vars)}"
      block.sequences.each do |seq|
        # TODO: i, q
        write_seq(seq)
      end
      f.puts
    end

    def write_seq(s)
      f.printf("%s %-20s %12d %2d %s %9d %s\n",
               s.empty? ? "e" : "s",
               s.source,
               s.start,
               s.size,
               s.strand,
               s.src_size,
               s.empty? ? s.status : s.text)
      if s.quality
        f.printf("q %-20s                           %s\n",
                 s.source, s.quality)
      end
      if s.i_data
        f.printf("i %-20s %s %s %s %s\n",
                 s.source, *s.i_data)
      end
    end
  end
  
end
