require 'spec_helper'

module Bio
  module MAF

    describe SQLiteIndex do

      describe ".build" do
        it "raises an error when trying to build an existing index" do
          expect {
            Bio::MAF::SQLiteIndex.build(nil, TestData + 'empty')
          }.to raise_error(/exists/)
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
            tuple[2].should == seq.start + seq.size
          end
        end

        describe "#index_tuple" do
          it "sets the bin correctly" do
            block = @p.parse_block
            seq = block.sequences.find { |s| s.source == @idx.sequence }
            tuple = @idx.index_tuple(block)
            # TODO: independently check this
            tuple[0].should == 1195
          end
        end
      end

    end
    
  end
  
end
