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

    end
    
  end
  
end
