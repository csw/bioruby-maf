require 'spec_helper'

module Bio::MAF

  describe Tiler do

    describe "#runs" do
      it "returns a uniform run properly" do
        a = Array.new(10, 'a')
        runs = Tiler.new.enum_for(:runs, a).to_a
        runs.should == [[0...10, 'a']]
      end
      it "yields a trailing item" do
        a = Array.new(10, 'a')
        a.fill('b', 8...10)
        runs = Tiler.new.enum_for(:runs, a).to_a
        runs.should == [[0...8, 'a'], [8...10, 'b']]
      end
      it "handles mixed contents" do
        spec = [[0...2, 'a'],
                [2...3, 'b'],
                [3...4, 'c'],
                [4...7, 'd']]
        a = Array.new(7, nil)
        spec.each { |range, obj| a.fill(obj, range) }
        runs = Tiler.new.enum_for(:runs, a).to_a
        runs.should == spec
      end
      it "handles overwrites" do
        spec = [[0...7, 'a'],
                [2...5, 'b'],
                [3...4, 'c'],
                [4...7, 'd']]
        a = Array.new(7, nil)
        spec.each { |range, obj| a.fill(obj, range) }
        runs = Tiler.new.enum_for(:runs, a).to_a
        runs.should == [[0...2, 'a'],
                        [2...3, 'b'],
                        [3...4, 'c'],
                        [4...7, 'd']]
      end
    end

  end

end
