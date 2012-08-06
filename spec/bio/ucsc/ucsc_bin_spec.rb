require 'spec_helper'

module Bio
  module Ucsc

    describe UcscBin do

      describe "#bin_from_range" do
        it "handles extended bin positions" do
          bin = UcscBin.bin_from_range(538457395, 538457395+44)
          bin.should == 13470
        end
      end

    end
    
  end
end
