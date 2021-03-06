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
      describe "#joinable_with?" do
        it "is false for blocks with different sequences" do
          p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          sp = %w(mm8 rn4 oryCun1 hg18 panTro2 rheMac2 canFam2 dasNov1 loxAfr1 echTel1)
          p.sequence_filter = { :only_species => sp }
          b1 = p.parse_block
          b2 = p.parse_block
          b1.joinable_with?(b2).should be_false
        end
        it "is true for blocks with same sequences" do
          p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          sp = %w(mm8 rn4 oryCun1 hg18 panTro2 rheMac2 canFam2 loxAfr1 echTel1)
          p.sequence_filter = { :only_species => sp }
          b1 = p.parse_block
          b2 = p.parse_block
          b1.joinable_with?(b2).should be_true
        end
      end
      describe "#to_bio_alignment" do
        it "returns a usable Bio::BioAlignment::Alignment" do
          p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          b = p.parse_block
          ba = b.to_bio_alignment
          ba.size.should == 10
          ba.sequences[0].id.should == "mm8.chr7"
          ba.sequences[0].seq.should =~ /^GGGCTGAGGGC--/
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

      describe "synteny data" do
        it "extracts basic data from i lines" do
          p = Parser.new(TestData + 'chr22_ieq2.maf',
                         :parse_extended => true)
          b = p.parse_block
          b.sequences[0].left_status_char.should be_nil
          b.sequences[0].left_status.should be_nil
          b.sequences[0].left_count.should be_nil
          b.sequences[0].right_status_char.should be_nil
          b.sequences[0].right_status.should be_nil
          b.sequences[0].right_count.should be_nil
          # works but let's not over-specify internal state
          #b.sequences[1].i_data.should == %w(N 0 C 0)
          b.sequences[1].left_status_char.should == 'N'
          b.sequences[1].left_status.should == :first
          b.sequences[1].right_status_char.should == 'C'
          b.sequences[1].right_status.should == :contiguous
          b.sequences[2].left_status.should == :contiguous
          b.sequences[2].right_status_char.should == 'I'
          b.sequences[2].right_status.should == :intervening
          b.sequences[2].right_count.should == 146
        end
      end

      describe "#to_bioalignment" do
        it "returns a usable Bio::BioAlignment::Sequence" do
          @parser = DummyParser.new
          line = "s human_unc 9077 8 + 10998 ACAGTATT"
          s = @parser.parse_seq_line(line, nil)
          as = s.to_bio_alignment
          as.id.should == "human_unc"
          as.seq.should == "ACAGTATT"
        end
      end

    end

  end
end
