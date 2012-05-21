require 'spec_helper'

module Bio
  module MAF

    describe Parser do

      describe "#open_file" do
        it "opens for reading the file passed as an argument"
        it "raises an error when the file does not exist"
        it "raises an error when the file is not in MAF format"
        it "parses the MAF header"
      end

      describe "#header" do
        it "raises an error when a MAF file is not open"
        it "provides version information"
        it "provides the scoring scheme"
        it "provides alignment parameters"
        it "provides arbitrary parameters"
      end

    end
    
  end
  
end
