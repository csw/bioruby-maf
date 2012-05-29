require 'bio-ucsc-api'

Given /^I have a region with start (\d+) and end (\d+)$/ do |r_start, r_end|
  @r_start = r_start.to_i
  @r_end = r_end.to_i
end

When /^I compute the smallest containing bin$/ do
  @bin = Bio::Ucsc::UcscBin.bin_from_range(@r_start, @r_end)
end

Then /^the bin should be (\d+)$/ do |expected_bin|
  @bin.should == expected_bin.to_i
end
