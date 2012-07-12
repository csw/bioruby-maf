require 'spec_helper'

module Bio::MAF

  describe Tiler do

    describe "#reference=" do
      it "takes a FASTARangeReader" do
        t = Tiler.new
        t.reference = FASTARangeReader.new((TestData + 'gap-sp1.fa').to_s)
        t.reference.read_interval(0, 3).should == 'CCA'
      end
      it "takes an open file" do
        (TestData + 'gap-sp1.fa').open do |fh|
          t = Tiler.new
          t.reference = fh
          t.reference.read_interval(0, 3).should == 'CCA'
        end
      end
      it "takes a FASTA literal" do
        ref = File.read(TestData + 'gap-sp1.fa')
        t = Tiler.new
        t.reference = ref
        t.reference.read_interval(0, 3).should == 'CCA'
      end
      it "takes a Pathname" do
        t = Tiler.new
        t.reference = TestData + 'gap-sp1.fa'
        t.reference.read_interval(0, 3).should == 'CCA'
      end
      it "takes a path string" do
        t = Tiler.new
        t.reference = (TestData + 'gap-sp1.fa').to_s
        t.reference.read_interval(0, 3).should == 'CCA'
      end
      it "takes a gzipped Pathname" do
        t = Tiler.new
        t.reference = TestData + 'gap-sp1.fa.gz'
        t.reference.read_interval(0, 3).should == 'CCA'
      end
      it "takes a gzipped path string" do
        t = Tiler.new
        t.reference = (TestData + 'gap-sp1.fa.gz').to_s
        t.reference.read_interval(0, 3).should == 'CCA'
      end
    end

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

  describe FASTARangeReader do
    describe "#read" do
      before(:each) do
        @r = FASTARangeReader.new('test/data/gap-sp1.fa')
        @s = 'CCAGGATGCTGGGCTGAGGGCAGTTGTGTCAGGGCGGTCCGGTGCAGGCA'
      end

      def check_range(z_start, z_end)
        @r.read_interval(z_start, z_end).should == @s.slice(z_start...z_end)
      end

      it "returns the entire sequence" do
        check_range(0, 50)
      end
      it "returns an entire line" do
        check_range(10, 20)
      end
      it "returns arbitrary components" do
        check_range(17, 41)
      end
    end
  end

end
      it "returns a sub-line fragment" do
        check_range(0, 4)
      end
