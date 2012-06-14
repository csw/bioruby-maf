When /^build an index on the reference sequence$/ do
  @idx = Bio::MAF::KyotoIndex.build(@parser, '%')
end

Given /^a Kyoto Cabinet index file "(.*?)"$/ do |name|
  @idx = Bio::MAF::KyotoIndex.open($test_data + name)
end

Then /^the index has at least (\d+) entries$/ do |size_spec|
  @idx.db.count.should be >= size_spec.to_i
end

When /^search for blocks between positions (\d+) and (\d+) of (\S+)$/ do |i_start, i_end, chr|
  int = Bio::GenomicInterval.zero_based(chr, i_start.to_i, i_end.to_i)
  @blocks = @idx.find([int], @parser, @block_filter).to_a
end

Then /^(\d+) blocks? (?:is|are) obtained$/ do |num|
  @blocks.size.should == num.to_i
end
