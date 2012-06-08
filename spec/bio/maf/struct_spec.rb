require 'spec_helper'

module Bio
  module MAF

    describe Struct do

      describe "#fmt" do
        it "presents all members in order" do
          @s = Struct.new([[:a, :uint16],
                           [:b, :uint32],
                           [:c, :uint32],
                           [:d, :uint8]])
          @s.fmt.should == "S>L>L>C"
        end
      end

      describe "#extractor_fmt" do
        it "handles uint8" do
          @s = Struct.new([[:marker, :uint8]])
          @s.extractor_fmt(:marker).should == "C"
        end
        it "handles uint16" do
          @s = Struct.new([[:a, :uint16]])
          @s.extractor_fmt(:a).should == "S>"
        end
        it "handles uint32" do
          @s = Struct.new([[:a, :uint32]])
          @s.extractor_fmt(:a).should == "L>"
        end
        it "handles uint64" do
          @s = Struct.new([[:a, :uint64]])
          @s.extractor_fmt(:a).should == "Q>"
        end
        it "skips uint8" do
          @s = Struct.new([[:dummy, :uint8],
                           [:a, :uint64]])
          @s.extractor_fmt(:a).should == "@1Q>"
        end
        it "skips uint16" do
          @s = Struct.new([[:dummy, :uint16],
                           [:a, :uint64]])
          @s.extractor_fmt(:a).should == "@2Q>"
        end
        it "skips uint32" do
          @s = Struct.new([[:dummy, :uint32],
                           [:a, :uint64]])
          @s.extractor_fmt(:a).should == "@4Q>"
        end
        it "skips uint64" do
          @s = Struct.new([[:dummy, :uint64],
                           [:a, :uint64]])
          @s.extractor_fmt(:a).should == "@8Q>"
        end
        it "extracts multiple leading elements" do
          @s = Struct.new([[:a, :uint16],
                           [:b, :uint32],
                           [:c, :uint32]])
          @s.extractor_fmt(:a, :b).should == "S>L>"
        end
        it "extracts multiple offset elements" do
          @s = Struct.new([[:a, :uint16],
                           [:b, :uint32],
                           [:c, :uint32],
                           [:d, :uint8]])
          @s.extractor_fmt(:b, :c).should == "@2L>L>"
        end
      end

      describe ""

    end

  end
end
