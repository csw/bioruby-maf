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

    describe Block do
      describe "#find_gaps" do
        it "finds a single 14-base gap" do
          p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          p.sequence_filter = { :only_species => %w(mm8 rn4 hg18 canFam2 loxAfr1) }
          block = p.parse_block
          gaps = block.find_gaps
          gaps.size.should == 1
          gaps[0][0].should == 34
          gaps[0][1].should == 14
        end
      end
      describe "#remove_gaps!" do
        it "removes a single 14-base gap" do
          p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          p.sequence_filter = { :only_species => %w(mm8 rn4 hg18 canFam2 loxAfr1) }
          block = p.parse_block
          block.sequences.size.should == 5
          block.text_size.should == 54
          block.remove_gaps!
          block.text_size.should == 40
        end
      end
    end

    describe Sequence do
      before(:each) do
        @parser = DummyParser.new
      end

      describe "#gapped?" do
        it "is false for sequences with no gaps" do
          line = "s human_unc 9077 8 + 10998 ACAGTATT"
          s = @parser.parse_seq_line(line, nil)
          s.gapped?.should be_false
        end
        it "is true for sequences with gaps" do
          line = "s human_unc 9077 8 + 10998 AC-AGTATT"
          s = @parser.parse_seq_line(line, nil)
          s.gapped?.should be_true
        end
      end

      describe "#text_range" do
        it "returns 0...text.size for a spanning interval" do
          line = "s human_unc 9077 8 + 10998 ACAGTATT"
          s = @parser.parse_seq_line(line, nil)
          range = s.text_range(9077...(9077 + 8))
          range.should == (0...(s.text.size))
        end
        it "returns 0...text.size for a gapped spanning interval" do
          line = "s human_unc 9077 8 + 10998 AC--AGTATT"
          s = @parser.parse_seq_line(line, nil)
          range = s.text_range(9077...(9077 + 8))
          range.should == (0...(s.text.size))
        end
        it "handles a leading subset" do
          line = "s human_unc 9077 8 + 10998 ACAGTATT"
          s = @parser.parse_seq_line(line, nil)
          range = s.text_range(9077...(9077 + 2))
          range.should == (0...2)
        end
        it "handles a trailing subset" do
          line = "s human_unc 9077 8 + 10998 ACAGTATT"
          s = @parser.parse_seq_line(line, nil)
          range = s.text_range(9079...9085)
          range.should == (2...8)
        end
        it "handles a gap in the middle" do
          line = "s human_unc 9077 8 + 10998 AC--AGTATT"
          s = @parser.parse_seq_line(line, nil)
          range = s.text_range(9078...(9077 + 8))
          range.should == (1...(s.text.size))
        end
        it "errors on a range starting before" do
          expect {
            line = "s human_unc 9077 8 + 10998 ACAGTATT"
            s = @parser.parse_seq_line(line, nil)
            range = s.text_range(9076...(9077 + 8))
          }.to raise_error
        end
        it "errors on a range ending after" do
          expect {
            line = "s human_unc 9077 8 + 10998 ACAGTATT"
            s = @parser.parse_seq_line(line, nil)
            range = s.text_range(9076...(9077 + 9))
          }.to raise_error
        end

      end

    end

  end
end
