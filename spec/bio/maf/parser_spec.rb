require 'spec_helper'

module Bio
  module MAF

    describe Header do
      before(:each) do
        @p = Parser.new(TestData + 't1.maf')
      end

      it "provides version information" do
        @p.header.version.should == '1'
      end
      it "provides the scoring scheme" do
        @p.header.scoring.should == 'humor.v4'
      end
      it "provides alignment parameters" do
        @p.header.alignment_params.should =~ /humor.v4 R=30/
      end

      it "presents multiline parameters correctly" do
        @p.header.alignment_params.should == "humor.v4 R=30 M=10 /cluster/data/hg15/bed/blastz.mm3/axtNet300/chr1.maf /cluster/data/hg15/bed/blastz.rn3/axtNet300/chr1.maf"
      end

      it "provides arbitrary parameters"
    end

    shared_examples "parsers" do
      
      describe "creation" do
        it "opens a file specified as a String argument"
        it "takes an IO object as an open file"
        it "raises an error when the file does not exist" do
          expect {
            described_class.new("/doesnotexist")
          }.to raise_error(Errno::ENOENT)
        end
        it "raises an error when the file is not in MAF format" do
          expect {
            described_class.new(TestData + '../../Rakefile')
          }.to raise_error
        end
      end

      describe "#header" do
        it "parses the MAF header" do
          p = described_class.new(TestData + 't1.maf')
          p.header.should_not be_nil
        end
      end

      context "at end of file" do
        describe "#parse_block" do
          it "returns nil"
        end
      end

      describe "#parse_block" do
        it "returns an alignment block" do
          p = described_class.new(TestData + 't1.maf')
          b = p.parse_block()
          b.should_not be_nil
        end
        it "raises an exception for malformed data"
      end

      it "handles absent alignment parameters" do
        p = described_class.new(TestData + 'chrY-1block.maf')
        b = p.parse_block()
        b.should_not be_nil
      end

      it "raises an exception on inconsistent sequence length" do
        pending
        ## can't just do string length, have to skip over hyphens
      end

    end

    describe Parser do
      include_examples "parsers"
    end

    describe ChunkParser do
      include_examples "parsers"

      def with_const_value(mod, sym, value)
        old = mod.const_get(sym)
        mod.const_set(sym, value)
        begin
          yield
        ensure
          mod.const_set(sym, old)
        end
      end

      it "sets last block position correctly" do
        p = ChunkParser.new(TestData + 'mm8_subset_a.maf')
        p.last_block_pos.should == 1103
      end

      it "parses larger files" do
        p = ChunkParser.new(TestData + 'mm8_chr7_tiny.maf')
        p.each_block { |block| block }
      end

      it "yields the correct number of blocks over chunk boundaries" do
        with_const_value(Bio::MAF::ChunkParser, :CHUNK_SIZE, 2048) do
          p = ChunkParser.new(TestData + 'mm8_chr7_tiny.maf')
          ref_scores = %w(10542.0 -33148.0 87527.0 185399.0 30120.0 58255.0 2607.0 8132.0)
          scores = []
          p.each_block do |block|
            scores << block.vars[:score]
          end
          scores.should == ref_scores
        end
      end

      it "sets last_block_pos correctly" do
        with_const_value(Bio::MAF::ChunkParser, :CHUNK_SIZE, 2048) do
          p = ChunkParser.new(TestData + 'mm8_chr7_tiny.maf')
          #p.parse_block
          p.last_block_pos.should == 1103
        end
      end

      it "handles sequence lines over chunk boundaries" do
        with_const_value(Bio::MAF::ChunkParser, :CHUNK_SIZE, 2048) do
          p = ChunkParser.new(TestData + 'mm8_chr7_tiny.maf')
          p.parse_block
          block = p.parse_block
          break_seq = block.raw_seq(4)
          break_seq.text.size.should == 156
        end
      end
    end

    describe LineReader do
      describe "creation" do
        it "takes a file-like object as argument" do
          expect {
            @f = File.open(TestData + 't1.maf', 'r')
            @reader = LineReader.new(@f)
            @f.close
          }.not_to raise_error
        end
      end

      context "with a file open" do
        before(:each) do
          @f = File.open(TestData + 't1.maf', 'r')
          @reader = LineReader.new(@f)
        end

        after(:each) do
          @f.close
        end

        describe "#next_line" do
          it "returns the next line of the file" do
            @reader.next_line.should =~ /^##maf/
          end
        end

        describe "#rewind" do
          it "causes #next_line to return the same line again" do
            l1 = @reader.next_line
            @reader.rewind
            @reader.next_line.should == l1
          end
        end
      end
    end
    
  end
  
end
