When /^build an index on the reference sequence in "(.*?)"$/ do |idx_path|
  pending("index")
  @idx_path = idx_path
  Bio::MAF::SQLiteIndex.build(@parser, idx_path)
end

Then /^the index has (\d+) entries$/ do |size_spec|
  idx = Bio::MAF::SQLiteIndex.open(idx_path)
  idx.size.should == size_spec.to_i
end
