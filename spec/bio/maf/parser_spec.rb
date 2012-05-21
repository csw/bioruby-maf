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

    describe Parser do

      describe "creation" do
        it "opens a file specified as a String argument"
        it "takes an IO object as an open file"
        it "raises an error when the filfe does not exist" do
          expect {
            Parser.new("/doesnotexist")
          }.to raise_error(Errno::ENOENT)
        end
        it "raises an error when the file is not in MAF format" do
          expect {
            Parser.new(TestData + '../../Rakefile')
          }.to raise_error
        end
      end

      describe "#header" do
        it "parses the MAF header" do
          p = Parser.new(TestData + 't1.maf')
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
          pending
          p = Parser.new(TestData + 't1.maf')
          b = p.parse_block()
          b.should_not be_nil
        end
        it "raises an exception for malformed data"
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
