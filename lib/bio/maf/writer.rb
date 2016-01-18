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
      f.write "##maf #{flatten_vars(header.vars)}\n"
      f.write "##{header.alignment_params}\n" if header.alignment_params
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
      lines << "\n"
      f.write(lines.join("\n"))
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

  FASTA_LINE_LEN = 72

  class FASTAWriter

    def initialize(outf)
      @f = outf
    end

    def write_block(block)
      block.sequences.each do |seq|
        write_sequence(seq) unless seq.empty?
      end
    end

    def write_sequence(seq)
      @f.puts(">#{seq.fasta_desc}")
      0.step(seq.text.size, FASTA_LINE_LEN) do |pos|
        @f.puts(seq.text.slice(pos, FASTA_LINE_LEN))
      end
    end

    def close
      @f.close
    end
  end
  
end
