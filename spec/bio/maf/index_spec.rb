require 'spec_helper'

module Bio
  module MAF

    describe KyotoIndex do

      describe ".build" do
        it "accepts '%' as a path for an in-memory DB" do
          expect {
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
            @p.f.close
            @idx.close
          }.not_to raise_error
        end
        it "accepts .kct paths"
        it "rejects other paths"
        context "mm8_chr7" do
          before(:each) do 
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
          end
          it "uses the first sequence appearing as the reference sequence" do
            @idx.index_sequences.to_a.should == [["mm8.chr7", 0]]
          end
          it "creates 8 index entries" do
            keys = @idx.db.match_prefix("\xFF\x00")
            keys.size.should == 8
          end
          it "stores the sequence IDs" do
            @idx.db.match_prefix("sequence:").size.should == 1
          end
          it "stores the sequence IDs" do
            @idx.db.get("sequence:mm8.chr7").should == "0"
          end
          after(:each) do
            @idx.db.close
          end
        end
      end

      describe ".open" do
 
      end

      describe "#find" do
        context "mm8_chr7" do
          before(:each) do
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
          end

          it "returns a block given a range contained in the block" do
            l = @idx.find([GenomicInterval.zero_based('mm8.chr7',
                                                      80082334,
                                                      80082338)],
                                @p)
            l.size.should == 1
            l[0].offset.should == 16
          end

          after(:each) do
            @idx.db.close
            @p.f.close
          end
        end
      end

      describe "#fetch_list" do
        context "mm8_chr7" do
          before(:each) do
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
          end
          it "returns a block spec given a range contained in the block" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082338)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns a block spec with correct size" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082338)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
            l[0][1].should == 1087 # block size
          end
          it "returns a block spec given its range exactly" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082368)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns specs for adjoining blocks given a range partially in each" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082360,
                                                            80082370)])
            l.size.should == 2
            l.collect { |e| e[0] }.should == [16, 1103]
          end
          it "returns a block spec given a range ending in it" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082330,
                                                            80082339)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns no block spec given a zero-based range ending at a block start" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082330,
                                                            80082334)])
            l.size.should == 0
          end
          it "returns a block spec given a range beginning in it" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083009,
                                                            80083220)])
            l.size.should == 1
            l[0][0].should == 10113 # block offset
          end
          it "returns no block spec given a range beginning at its end" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083156,
                                                            80083200)])
            l.size.should == 0
          end
          it "returns specs for all blocks given a range fitting a larger bin" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            0,
                                                            80083200)])
            l.size.should == 8
          end
          it "returns no blocks given a range outside" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083200,
                                                            80083300)])
          end
          after(:each) do
            if @idx
              @idx.db.close
            end
          end
        end
      end

      describe "#entries_for" do
        before(:each) do
          @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          @block = @p.parse_block
          @idx = KyotoIndex.new('%')
        end
        context "single ref seq" do
          before(:each) do
            @idx.index_sequences = { 'mm8.chr7' => 0 }
            @e = @idx.entries_for(@block)
          end
          it "returns a two-element array" do
            @e[0].size.should == 2
          end
          it "gives the correct key data" do
            _, seq, bin, i_start, i_end = @e[0][0].unpack("CCS>L>L>")
            seq.should == 0
            bin.should == 1195
            i_start.should == 80082334
            i_end.should == 80082368
          end
          it "gives the correct offset" do
            b_offset, b_len = @e[0][1].unpack("Q>Q>")
            b_offset.should == 16
          end
          it "gives the correct length" do
            b_offset, b_len = @e[0][1].unpack("Q>Q>")
            b_len.should == 1087
          end
        end
        after(:each) do
          @p.f.close
          @idx.db.close
        end
      end

    end

  end # module MAF
  
end # module Bio
