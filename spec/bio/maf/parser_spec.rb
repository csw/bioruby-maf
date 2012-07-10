require 'spec_helper'

module Bio
  module MAF

    describe ParseContext do
      it "tracks the last block position"
    end

    describe ChunkReader do
      before(:each) do
        @f = (TestData + 'mm8_chr7_tiny.maf').open
      end
      describe "#initialize" do
        it "rejects a chunk size of zero" do
          expect {
            ChunkReader.new(@f, 0)
          }.to raise_error(/Invalid chunk size/)
        end
        it "rejects a negative chunk size" do
          expect {
            ChunkReader.new(@f, 0)
          }.to raise_error(/Invalid chunk size/)
        end
        it "rejects a chunk size not a power of 2" do
          expect {
            ChunkReader.new(@f, 1000)
          }.to raise_error(/Invalid chunk size/)
        end
        it "accepts a 4k chunk size" do
          expect {
            ChunkReader.new(@f, 4096)
          }.not_to raise_error
        end
        it "accepts an 8M chunk size" do
          expect {
            ChunkReader.new(@f, 8 * 1024 * 1024)
          }.not_to raise_error
        end
      end
      context "with 1K ChunkReader" do
        before(:each) do
          @r = ChunkReader.new(@f, 1024)
        end
 
        describe "#chunk_size=" do
          it "sets the chunk size" do
            @r.chunk_size = 8192
            @r.chunk_size.should == 8192
          end
          # it "sets the chunk shift" do
          #   @r.chunk_size = 8192
          #   @r.chunk_shift.should == 13
          # end
        end

        describe "#read_chunk" do
          it "returns a chunk of the specified length" do
            @r.read_chunk.bytesize == 1024
          end
          it "starts at position 0" do
            @r.pos.should == 0
          end
          it "advances the position" do
            @r.read_chunk
            @r.pos.should == 1024
          end
        end

        describe "#read_chunk_at" do
          it "returns data starting at the specified offset" do
            c = @r.read_chunk_at(59)
            c.start_with?("80082334").should be_true
          end
          it "handles a read starting exactly at a chunk boundary" do
            c = @r.read_chunk_at(1024)
            c.start_with?("   594").should be_true
          end
        end
      end
      after(:each) do
        @f.close
      end
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

      describe "#fetch_blocks" do
        shared_examples_for "any chunk size" do
          it "parses a single block" do
            fl = [[16, 1087]]
            blocks = @p.fetch_blocks(fl).to_a
            blocks.size.should == 1
            blocks[0].offset.should == 16
          end
          it "parses several consecutive blocks" do
            fl = [[16, 1087], [1103, 1908], [3011, 2027]]
            blocks = @p.fetch_blocks(fl).to_a
            blocks.size.should == 3
            blocks.collect {|b| b.offset}.should == [16, 1103, 3011]
          end
          it "parses consecutive blocks further ahead" do
            fl = [[5038, 1647], [6685, 829]]
            blocks = @p.fetch_blocks(fl).to_a
            blocks.size.should == 2
            blocks.collect {|b| b.offset}.should == [5038, 6685]
          end
          it "parses nonconsecutive blocks" do
            fl = [[16, 1087], [3011, 2027]]
            blocks = @p.fetch_blocks(fl).to_a
            blocks.size.should == 2
            blocks.collect {|b| b.offset}.should == [16, 3011]
          end
          it "takes a block argument" do
            fl = [[16, 1087], [1103, 1908], [3011, 2027]]
            n = 0
            @p.fetch_blocks(fl) do |block|
              n += 1
            end
            n.should == 3
          end
        end
        context "with 4K chunk size" do
          before(:each) do
            @p = described_class.new(TestData + 'mm8_chr7_tiny.maf',
                                     :chunk_size => 4096,
                                     :random_chunk_size => 4096)
          end
          it_behaves_like "any chunk size"
        end
        context "with 1K chunk size" do
          before(:each) do
            @p = described_class.new(TestData + 'mm8_chr7_tiny.maf',
                                     :chunk_size => 1024,
                                     :random_chunk_size => 1024)
          end
          it_behaves_like "any chunk size"
        end
        context "after parsing to end" do
          before(:each) do
            @p = described_class.new(TestData + 'mm8_chr7_tiny.maf',
                                     :chunk_size => 4096,
                                     :random_chunk_size => 4096)
            @p.each_block { |b| nil }
          end
          it_behaves_like "any chunk size"
        end
        context "with 8M chunk size" do
          before(:each) do
            @p = described_class.new(TestData + 'mm8_chr7_tiny.maf',
                                     :chunk_size => 8 * 1024 * 1024,
                                     :random_chunk_size => 8 * 1024 * 1024)
          end
          it_behaves_like "any chunk size"
        end
        after(:each) do
          @p.f.close
        end
      end

      describe "#each_block" do
        it "returns an Enumerator when called without a block" do
          p = described_class.new(TestData + 'mm8_chr7_tiny.maf')
          p.each_block.count.should == 8
        end
      end

      describe "sequence_filter" do
        before(:each) do
          @p = described_class.new(TestData + 'mm8_mod_a.maf')
        end
        it "restricts sequences parsed" do
          @p.sequence_filter = { :only_species => %w(mm8 rn4) }
          @p.parse_block.sequences.size.should == 2
        end
        it "matches at the species delimiter rather than a prefix" do
          @p.sequence_filter = { :only_species => %w(mm8 hg18) }
          @p.parse_block.sequences.size.should == 2
        end
        it "sets filtered? when modified" do
          @p.sequence_filter = { :only_species => %w(mm8 rn4) }
          @p.parse_block.filtered?.should be_true
        end
        it "does not set filtered? when unmodified" do
          @p.sequence_filter = {
            :only_species => %w(mm8 rn4 oryCun1 hg18 hg181)
          }
          @p.parse_block.filtered?.should be_false
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

      it "gives the correct number of sequences" do
        p = described_class.new(TestData + 'mm8_chr7_tiny.maf')
        block = p.parse_block
        block.sequences.size.should == 10
      end

      it "handles absent alignment parameters" do
        p = described_class.new(TestData + 'chrY-1block.maf')
        b = p.parse_block()
        b.should_not be_nil
      end

      it "parses larger files" do
        p = described_class.new(TestData + 'mm8_chr7_tiny.maf')
        expect {
          p.each_block { |block| block }
        }.not_to raise_error
      end

      it "handles trailing comments" do
        p = described_class.new(TestData + 't1a.maf')
        expect {
          p.each_block { |block| block }
        }.not_to raise_error
      end

      it "raises an exception on inconsistent sequence length" do
        pending
        ## can't just do string length, have to skip over hyphens
      end

      it "tracks block start offsets correctly" do
        pa = []
        p = described_class.new(TestData + 'mm8_chr7_tiny.maf')
        p.each_block { |b| pa << b.offset }
        pa.should == [16, 1103, 3011, 5038, 6685, 7514, 9022, 10113]
      end

      it "reports block sizes correctly" do
        p = described_class.new(TestData + 'mm8_chr7_tiny.maf')
        block = p.parse_block
        block.size.should == 1087
      end

      it "parses very large blocks" do
        p = described_class.new(TestData + 'big-block.maf')
        n = 0
        p.each_block { |b| n += 1 }
        n.should == 490
      end

    end

    describe Parser do
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

      describe "#merge_fetch_list" do
        before(:each) do
          @p = described_class.new(TestData + 'mm8_chr7_tiny.maf')
        end
        it "passes through single records" do
          fl = [[16, 1087]]
          @p.merge_fetch_list(fl).should == [[16, 1087, [16]]]
        end
        it "passes through non-contiguous records" do
          fl = [[16, 1087], [3011, 2027]]
          @p.merge_fetch_list(fl).should == [[16, 1087, [16]],
                                             [3011, 2027, [3011]]]
        end
        it "merges contiguous records" do
          fl = [[16, 1087], [1103, 1908], [3011, 2027]]
          @p.merge_fetch_list(fl).should == [[16, 5022, [16, 1103, 3011]]]
        end
        after(:each) do
          @p.f.close
        end
      end

      it "sets last block position correctly" do
        p = Parser.new(TestData + 'mm8_subset_a.maf')
        p.last_block_pos.should == 1103
      end

      context "with 2k chunk size" do
        before(:each) do
          @p = Parser.new(TestData + 'mm8_chr7_tiny.maf',
                          :chunk_size => 2048)
        end
        it "yields the correct number of blocks over chunk boundaries" do
          ref_scores = %w(10542.0 -33148.0 87527.0 185399.0
                          30120.0 58255.0 2607.0 8132.0)
          scores = []
          @p.each_block do |block|
            scores << block.vars[:score]
          end
          scores.should == ref_scores
        end
        it "sets last_block_pos correctly" do
          @p.last_block_pos.should == 1103
        end
        it "handles sequence lines over chunk boundaries" do
          @p.parse_block
          block = @p.parse_block
          break_seq = block.raw_seq(4)
          break_seq.text.size.should == 156
        end

        it "tracks block start offsets correctly over chunk boundaries" do
          pa = []
          @p.each_block { |b| pa << b.offset }
          pa.should == [16, 1103, 3011, 5038, 6685, 7514, 9022, 10113]
        end
      end

    end

  end
  
end
