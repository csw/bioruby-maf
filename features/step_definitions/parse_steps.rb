When /^I open it with a MAF reader$/ do
  @opts ||= {}
  @parser = Bio::MAF::Parser.new(@src_f, @opts)
end

When /^I enable the :(\S+) parser option$/ do |opt_s|
  if @parser
    opts = @parser.opts
  elsif @access
    opts = @access.parse_options
  else
    @opts ||= {}
    opts = @opts
  end
  opts[opt_s.to_sym] = true
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
  @block.sequences.size.should == n_seq.to_i
end

Then /^block (\d+) has (\d+) sequences$/ do |block_n, n_seq|
  @blocks[block_n.to_i].sequences.size.should == n_seq.to_i
end

Then /^sequence (\d+) has (\w.*?) "(.*?)"$/ do |i, method, str|
  method_sym = method.gsub(/ /, '_').to_sym
  @block.raw_seq(i.to_i).send(method_sym).should == str
end

Then /^sequence (\d+) has (\w.*?) (\d+)\s*$/ do |i, method, num|
  method_sym = method.gsub(/ /, '_').to_sym
  @block.raw_seq(i.to_i).send(method_sym).should == num.to_i
end

Then /^sequence (\d+) has (\w.*?) :(\S+)\s*$/ do |i, method, sym_s|
  method_sym = method.gsub(/ /, '_').to_sym
  value_sym = sym_s.to_sym
  @block.raw_seq(i.to_i).send(method_sym).should == value_sym
end

Then /^sequence (\S+) of block (\d+) has (\w.*?) "(.*?)"$/ do |chr, i, method, str|
  seq = @blocks[i.to_i].sequences.find { |seq| seq.source == chr }
  method_sym = method.gsub(/ /, '_').to_sym
  seq.send(method_sym).should == str
end

Then /^sequence (\S+) of block (\d+) has (\w.*?) (\d+)$/ do |chr, i, method, num|
  seq = @blocks[i.to_i].sequences.find { |seq| seq.source == chr }
  method_sym = method.gsub(/ /, '_').to_sym
  seq.send(method_sym).should == num.to_i
end

Then /^sequence (\S+) of block (\d+) has (\w.*?) :(\S+)$/ do |chr, i, method, sym_s|
  seq = @blocks[i.to_i].sequences.find { |seq| seq.source == chr }
  method_sym = method.gsub(/ /, '_').to_sym
  seq.send(method_sym).should == sym_s.to_sym
end
