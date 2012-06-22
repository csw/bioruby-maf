module Bio::MAF

  class Writer
    attr_reader :f, :path

    def initialize(fspec)
      if fspec.respond_to? :write
        @f = fspec
        if fspec.respond_to? :path
          @path = fspec.path
        end
      else
        @path = fspec
        @f = File.open(fspec, 'w')
      end
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
      lines = ["a #{flatten_vars(block.vars)}"]
      block.sequences.each do |seq| 
        write_seq(seq, lines)
      end
      lines << " "
      f.puts lines.join("\n")
    end

    def write_seq(s, lines)
      lines << sprintf("%s %-20s %12d %2d %s %9d %s",
                       s.empty? ? "e" : "s",
                       s.source,
                       s.start,
                       s.size,
                       s.strand,
                       s.src_size,
                       s.empty? ? s.status : s.text)
      if s.quality
        lines << sprintf("q %-20s                           %s",
                         s.source, s.quality)
      end
      if s.i_data
        lines << sprintf("i %-20s %s %s %s %s",
                         s.source, *s.i_data)
      end
    end
  end
  
end
