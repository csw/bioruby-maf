When /^slice the resulting blocks according to the given interval$/ do
  # @blocks and @interval
  @blocks = @blocks.collect { |b| b.slice(@interval) }
end

When /^I extract a slice over the genomic interval$/ do |table|
  # table is a Cucumber::Ast::Table
  intervals = table.hashes.collect do |row|
    Bio::GenomicInterval.zero_based(row['chrom'],
                                    row['start'].to_i,
                                    row['end'].to_i)
  end
  intervals.size.should == 1
  @blocks = @access.slice(intervals[0])
end
