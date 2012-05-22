Given /^in a temp file$/ do
  @src_f = Tempfile.new(['rspec', '.maf'])
  @src_f.write(@src)
  @src_f.close
end

When /^I open it with a MAF reader$/ do
  @parser = Bio::MAF::Parser.new(@src_f)
end

When /^I open it with a chunked MAF reader$/ do
  @parser = Bio::MAF::ChunkParser.new(@src_f)
end

Then /^the MAF version should be "(.*?)"$/ do |v_spec|
  @parser.header.version.to_s.should == v_spec
end

Then /^the scoring scheme should be "(.*?)"$/ do |s_spec|
  @parser.header.scoring.should == s_spec
end

Then /^the alignment parameters should be "(.*?)"$/ do |a_spec|
  @parser.header.alignment_params.should == a_spec
end

Then /^an alignment block can be obtained$/ do
  @block = @parser.parse_block
  @block.should_not be_nil
end

Then /^the alignment block has (\d+) sequences$/ do |n_seq|
  @block.size.should == n_seq.to_i
end

Then /^sequence (\d+) has source "(.*?)"$/ do |i, src|
  @block.raw_seq(i.to_i).source.should == src
end

Then /^sequence (\d+) has start (\d+)$/ do |i, start|
  @block.raw_seq(i.to_i).start.should == start.to_i
end

Then /^sequence (\d+) has size (\d+)$/ do |i, size|
  @block.raw_seq(i.to_i).size.should == size.to_i
end

Then /^sequence (\d+) has strand "(.*?)"$/ do |i, strand|
  @block.raw_seq(i.to_i).strand.should == strand.to_sym
end

Then /^sequence (\d+) has source size (\d+)$/ do |i, src_size|
  @block.raw_seq(i.to_i).src_size.should == src_size.to_i
end

Then /^sequence (\d+) has text "(.*?)"$/ do |i, text|
  @block.raw_seq(i.to_i).text.should == text
end
