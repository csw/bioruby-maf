require 'spec_helper'

module Bio
  module MAF

    describe KyotoIndex do

      describe ".build" do
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
          after(:each) do
            @idx.db.close
          end
        end
      end

      describe "#fetch_list" do
        context "mm8_chr7" do
          before(:each) do
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
          end
          it "returns a block given a range contained in the block" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082338)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns a block given its range exactly" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082368)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns adjoining blocks given a range partially in each" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082360,
                                                            80082370)])
            l.size.should == 2
            l.collect { |e| e[0] }.should == [16, 1103]
          end
          it "returns a block given a range ending in it" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082330,
                                                            80082339)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns no block given a zero-based range ending at its start" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082330,
                                                            80082334)])
            l.size.should == 0
          end
          it "returns a block given a range beginning in it" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083009,
                                                            80083220)])
            l.size.should == 1
            l[0][0].should == 10113 # block offset
          end
          it "returns no block given a range beginning in at its end" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083156,
                                                            80083200)])
            l.size.should == 0
          end
          it "returns all blocks given a range fitting a larger bin" do
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

    describe SQLiteIndex do

      describe ".open" do
        it "raises an error on a nonexistent index" do
          expect {
            @idx = SQLiteIndex.open(TestData + 'no_such_file')
          }.to raise_error(/does not exist/)
        end
        it "raises an error on an index without name metadata" do
          expect {
            @idx = SQLiteIndex.open(TestData + 'empty.db')
          }.to raise_error(/not a usable index database/)
        end
        it "raises an error on an index DB without ref seq data" do
          expect {
            @idx = SQLiteIndex.open(TestData + 'mm8_chr7_tiny.nometa.index')
          }.to raise_error(/not a usable index database/)
        end
        it "sets the sequence name correctly" do
          @idx = SQLiteIndex.open(TestData + 'mm8_chr7_tiny.index')
          @idx.sequence.should == "mm8.chr7"
        end
        it "sets the table name correctly" do
          @idx = SQLiteIndex.open(TestData + 'mm8_chr7_tiny.index')
          @idx.table_name.should == "seq_mm8_chr7"
        end
        after(:each) do
          if @idx
            @idx.db.disconnect
          end
        end
      end

      describe "#fetch_list" do
        context "mm8_chr7" do
          before(:each) do
            @idx = SQLiteIndex.open(TestData + 'mm8_chr7_tiny.index')
          end
          it "returns a block given a range contained in the block" do
            l = @idx.fetch_list([80082334..80082338])
            l.size.should == 1
            l[0][0].should == 16
          end
          it "returns a block given its range exactly"
          it "returns adjoining blocks given a range partially in each"
          it "returns a block given a range ending in it"
          it "returns a block given a range beginning in it"
          it "returns no blocks given a range outside"
          after(:each) do
            if @idx
              @idx.db.disconnect
            end
          end
        end
      end

      describe ".build" do
        it "raises an error when trying to build an existing index" do
          expect {
            Bio::MAF::SQLiteIndex.build(nil, TestData + 'empty')
          }.to raise_error(/exists/)
        end
 
        context "mm8_chr7" do
          before(:each) do 
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = Bio::MAF::SQLiteIndex.build(@p, ":memory:")
          end
          after(:each) do
            @idx.db.disconnect
          end

          it "uses the first sequence appearing as the reference sequence" do
            @idx.sequence.should == "mm8.chr7"
          end
          it "creates an index table named after the ref sequence" do
            @idx.count_tables("seq_mm8_chr7").should == 1
          end
          it "creates 8 index rows" do
            row = @idx.db.select_one("select count(*) from seq_mm8_chr7")
            count = row[0].to_i
            count.should == 8
          end
          it "creates an index" do
            r = @idx.db.select_one(<<-EOF)
SELECT COUNT(*) FROM sqlite_master
WHERE name = 'idx_seq_mm8_chr7'
AND type = 'index'
EOF
            r[0].to_i.should == 1
          end
        end
      end

      describe "#count_tables" do
        before(:each) do 
          @idx = SQLiteIndex.new(":memory:")
        end
        after(:each) do
          @idx.db.disconnect
        end
        it "reports no tables if none exist" do
          @idx.count_tables("xyzzy").should == 0
        end
        it "reports one table correctly" do
          @idx.db.do("create table abc (a numeric(4) primary key)")
          @idx.count_tables("abc").should == 1
        end
        it "does LIKE queries given a % escape" do
          @idx.db.do("create table abc1 (a numeric(4) primary key)")
          @idx.db.do("create table abc2 (a numeric(4) primary key)")
          @idx.count_tables("abc%").should == 2
        end
      end

      describe "#create_schema" do
        it "creates a metadata table" do
          idx = SQLiteIndex.new(":memory:")
          idx.create_schema
          idx.count_tables("metadata").should == 1
        end
      end

      describe "with mm8_chr7 data" do
        before(:each) do 
          @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          @idx = SQLiteIndex.new(":memory:")
          @idx.sequence = "mm8.chr7"
        end
        after(:each) do
          @idx.db.disconnect
        end

        describe "#index_tuple" do
          it "preserves the start and end" do
            block = @p.parse_block
            seq = block.sequences.find { |s| s.source == @idx.sequence }
            tuple = @idx.index_tuple(block)
            tuple[1].should == seq.start
            tuple[2].should == seq.start + seq.size - 1
          end
          it "sets the bin correctly" do
            block = @p.parse_block
            seq = block.sequences.find { |s| s.source == @idx.sequence }
            tuple = @idx.index_tuple(block)
            # TODO: independently check this
            tuple[0].should == 1195
          end
          it "includes the block position" do
            block = @p.parse_block
            seq = block.sequences.find { |s| s.source == @idx.sequence }
            tuple = @idx.index_tuple(block)
            # TODO: independently check this
            tuple[3].should == 16
          end
        end
      end

    end # SQLiteIndex
    
  end # module MAF
  
end # module Bio
