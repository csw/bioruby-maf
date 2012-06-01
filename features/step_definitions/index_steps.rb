When /^build an index on the reference sequence$/ do
  @idx = Bio::MAF::KyotoIndex.build(@parser, '%')
end

Then /^the index has at least (\d+) entries$/ do |size_spec|
  @idx.db.count.should be >= size_spec.to_i
end

When /^search for blocks between positions (\d+) and (\d+) of (\S+)$/ do |i_start, i_end, chr|
  int = Bio::GenomicInterval.zero_based(chr, i_start.to_i, i_end.to_i)
  @blocks = @idx.find([int], @parser)
end

Then /^(\d+) blocks are obtained$/ do |num|
  @blocks.size.should == num.to_i
end
