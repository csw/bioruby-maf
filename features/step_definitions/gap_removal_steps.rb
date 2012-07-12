Then /^the alignment block is marked as filtered$/ do
  @block.filtered?.should be_true
end

Then /^(\d+) gaps? (?:is|are) found with length \[(\d+)\]$/ do |n_gaps, gap_sizes_s|
  gaps = @block.find_gaps
  gaps.size.should == n_gaps.to_i
  e_gap_sizes = gap_sizes_s.split(/,\s*/).collect { |n| n.to_i }
  gap_sizes = gaps.collect { |gap| gap[1] }
  gap_sizes.should == e_gap_sizes
end

When /^gaps are removed$/ do
  @block.remove_gaps!
end

Then /^the text size of the block is (\d+)$/ do |e_text_size|
  @block.text_size.should == e_text_size.to_i
end

Then /^the text size of block (\d+) is (\d+)$/ do |n, e_text_size|
  @blocks[n.to_i].text_size.should == e_text_size.to_i
end
