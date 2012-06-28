require 'spec_helper'

module Bio
  module MAF

    describe KyotoIndex do
      def has_at_least_n_with_prefix(n, start)
        @idx.db.cursor_process do |cur|
          i = 0
          cur.jump(start)
          k = cur.get_key(true)
          while k && k.start_with?(start) && i < n
            i += 1
          end
          return i == n
        end
      end

      describe ".build" do
        it "accepts '%' as a path for an in-memory DB" do
          expect {
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
            @p.f.close
            @idx.close
          }.not_to raise_error
        end
        it "accepts .kct paths"
        it "rejects other paths"
        context "mm8_chr7" do
          before(:each) do 
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
          end
          it "uses the first sequence appearing as the reference sequence" do
            @idx.index_sequences.to_a.should == [["mm8.chr7", 0]]
          end
          it "creates 8 index entries" do
            has_at_least_n_with_prefix(8, "\xFF\x00").should be_true
          end
          it "stores the sequence IDs" do
            @idx.db.match_prefix("sequence:").size.should == 1
          end
          it "stores the sequence IDs" do
            @idx.db.get("sequence:mm8.chr7").should == "0"
          end
          describe "loads sequence data correctly" do
            before(:each) { @idx = @idx.reopen }
            it "uses the first sequence appearing as the reference sequence" do
              @idx.index_sequences.to_a.should == [["mm8.chr7", 0]]
            end
          end
          after(:each) do
            @idx.db.close
          end
        end
      end

      describe ".open" do
        it "opens an existing index successfully" do
          @idx = KyotoIndex.open(TestData + 'mm8_chr7_tiny.kct')
          @idx.db.count.should be > 8
        end
        it "populates #index_sequences" do
          @idx = KyotoIndex.open(TestData + 'mm8_chr7_tiny.kct')
          @idx.index_sequences.size.should be > 0
          @idx.index_sequences['mm8.chr7'].should == 0
        end
        after(:each) do
          @idx.db.close if @idx
        end
      end

      describe "#find" do
        context "mm8_chr7" do
          before(:each) do
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
          end

          it "returns a block given a range contained in the block" do
            l = @idx.find([GenomicInterval.zero_based('mm8.chr7',
                                                      80082334,
                                                      80082338)],
                                @p).to_a
            l.size.should == 1
            l[0].offset.should == 16
          end

          after(:each) do
            @idx.db.close
            @p.f.close
          end
        end
      end

      describe "#fetch_list" do
        context "mm8_chr7" do
          before(:each) do
            @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
            @idx = KyotoIndex.build(@p, '%')
          end
          it "returns a block spec given a range contained in the block" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082338)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns a block spec with correct size" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082338)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
            l[0][1].should == 1087 # block size
          end
          it "returns a block spec given its range exactly" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082334,
                                                            80082368)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns specs for adjoining blocks given a range partially in each" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082360,
                                                            80082370)])
            l.size.should == 2
            l.collect { |e| e[0] }.should == [16, 1103]
          end
          it "returns a block spec given a range ending in it" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082330,
                                                            80082339)])
            l.size.should == 1
            l[0][0].should == 16 # block offset
          end
          it "returns no block spec given a zero-based range ending at a block start" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80082330,
                                                            80082334)])
            l.size.should == 0
          end
          it "returns a block spec given a range beginning in it" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083009,
                                                            80083220)])
            l.size.should == 1
            l[0][0].should == 10113 # block offset
          end
          it "returns no block spec given a range beginning at its end" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083156,
                                                            80083200)])
            l.size.should == 0
          end
          it "returns specs for all blocks given a range fitting a larger bin" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            0,
                                                            80083200)])
            l.size.should == 8
          end
          it "returns no blocks given a range outside" do
            l = @idx.fetch_list([GenomicInterval.zero_based('mm8.chr7',
                                                            80083200,
                                                            80083300)])
          end
          after(:each) do
            if @idx
              @idx.db.close
            end
          end
        end
      end

      describe "#overlaps?" do
        before(:each) do
          @idx = KyotoIndex.new('%')
        end
        def check_overlap(x, y)
          i = x[0]...x[1]
          @idx.overlaps?(i, y[0], y[1])
        end
        it "handles equal intervals" do
          check_overlap([0, 10],
                        [0, 10]).should be_true
        end
        it "handles X contains Y" do
          check_overlap([0, 10],
                        [0, 9]).should be_true
          check_overlap([0, 10],
                        [1, 9]).should be_true
          check_overlap([0, 10],
                        [1, 10]).should be_true
        end
        it "handles Y contains X" do
          check_overlap([0, 9],
                        [0, 10]).should be_true
          check_overlap([1, 9],
                        [0, 10]).should be_true
          check_overlap([1, 10],
                        [0, 10]).should be_true
        end
        it "handles partial overlap" do
          check_overlap([0, 9],
                        [1, 10]).should be_true
          check_overlap([1, 10],
                        [0, 9]).should be_true
        end
        it "handles end cases" do
          check_overlap([0, 10],
                        [10, 15]).should be_false
          check_overlap([10, 15],
                        [0, 10]).should be_false
        end
        it "handles separated intervals" do
          check_overlap([0, 10], [15, 20]).should be_false
          check_overlap([15, 20], [0, 10]).should be_false
        end
        after(:each) do
          @idx.db.close
        end
      end

      describe "#entries_for" do
        before(:each) do
          @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
          @block = @p.parse_block
          @idx = KyotoIndex.new('%')
        end
        context "single ref seq" do
          before(:each) do
            @idx.index_sequences = { 'mm8.chr7' => 0 }
            @e = @idx.entries_for(@block)
          end
          it "gives the correct key data" do
            _, seq, bin, i_start, i_end = @e.keys.first.unpack("CCS>L>L>")
            seq.should == 0
            bin.should == 1195
            i_start.should == 80082334
            i_end.should == 80082368
          end
          it "gives the correct offset" do
            b_offset, b_len = @e.values.first.unpack("Q>L>")
            b_offset.should == 16
          end
          it "gives the correct length" do
            b_offset, b_len = @e.values.first.unpack("Q>L>")
            b_len.should == 1087
          end
        end
        after(:each) do
          @p.f.close
          @idx.db.close
        end
      end

    end

    describe "#species" do
      before(:each) do
        @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
        @idx = KyotoIndex.build(@p, '%')
      end
      shared_examples "species" do
        it "records the correct number of species" do
          @idx.species.size.should == 11
        end
        it "sets species_max_id correctly" do
          @idx.species_max_id.should == 10
        end
      end
      describe "after building index" do
        include_examples "species"
        it "records species in order" do
          @idx.db["species:mm8"].should == "0"
        end
      end
      describe "after loading index" do
        before(:each) { @idx = @idx.reopen }
        include_examples "species"
      end
    end

    describe "Filter classes" do
      before(:each) do 
        @p = Parser.new(TestData + 'mm8_chr7_tiny.maf')
        @idx = KyotoIndex.build(@p, '%')
      end

      describe AllSpeciesFilter do
        def fake_entry_with(species_l)
          ids = species_l.collect {|s| @idx.species.fetch(s)}
          vec = ids.collect { |id| 1 << id }.reduce(0, :|)
          return ['', [0, 0, 0, 0, vec].pack(KyotoIndex::VAL_FMT)]
        end

       context "with an empty set" do
          before(:each) do
            @filter = AllSpeciesFilter.new([], @idx)
          end
          it "matches anything" do
            e = fake_entry_with(%w(mm8 rn4 oryCun1))
            @filter.match(e).should be_true
          end
        end
        context "with [mm8 rn4]" do
          before(:each) do
            @filter = AllSpeciesFilter.new(%w(mm8 rn4), @idx)
          end
          it "does not match an empty entry" do
            e = fake_entry_with(%w())
            KVHelpers.extract_species_vec(e).should == 0
            @filter.bs.should_not == 0
            @filter.match(e).should be_false
          end
          it "does not match an entry with mm8" do
            e = fake_entry_with(%w(mm8))
            @filter.match(e).should be_false
          end
          it "does not match an entry with mm8 oryCun1" do
            e = fake_entry_with(%w(mm8 oryCun1))
            @filter.match(e).should be_false
          end
          it "matches an entry with mm8 rn4" do
            e = fake_entry_with(%w(mm8 rn4))
            @filter.match(e).should be_true
          end
          it "does not match an entry with mm8 rn4 oryCun1" do
            e = fake_entry_with(%w(mm8 rn4 oryCun1))
            @filter.match(e).should be_true
          end
        end
      end # AllSpeciesFilter

      describe AtLeastNSequencesFilter do
        def fake_entry_with(n)
          return ['', [0, 0, 0, n, 0].pack(KyotoIndex::VAL_FMT)]
        end
        context "n = 3" do
          before(:each) do
            @filter = AtLeastNSequencesFilter.new(3, @idx)
          end
          it "does not match 2 sequences" do
            e = fake_entry_with(2)
            @filter.match(e).should be_false
          end
          it "matches 3 sequences" do
            e = fake_entry_with(3)
            @filter.match(e).should be_true
          end
        end
      end # AtLeastNSequencesFilter
      
      after(:each) do
        @idx.close
      end
    end # filter classes

  end # module MAF
  
end # module Bio
