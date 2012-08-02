Given /^a MAF source file "(.*?)"$/ do |src|
  @src_f = $test_data + src
  @src_f.exist?.should be_true
end

Given /^MAF data:$/ do |string|
  @src_f = Tempfile.new(['rspec', '.maf'])
  @src_f.write(string)
  @src_f.rewind
end

When /^I select FASTA output$/ do
  @dst = Tempfile.new(['cuke', ".#{@out_fmt.to_s}"])
  @writer = Bio::MAF::FASTAWriter.new(@dst)
end

When /^process the file$/ do
  @parser.each_block do |block|
    @writer.write_block(block)
  end
  @writer.close
end

Then /^the output should match "(.*?)"$/ do |ref|
  ref_p = $test_data + ref
  ref_p.exist?.should be_true
  #system("diff #{ref} #{@dst.path} >/dev/null 2>&1").should be_true
  File.read(@dst.path).should == File.read(ref_p)
end

Then /^the output should be:$/ do |string|
  File.read(@dst.path).should == string
end

After do
  if @dst
    @dst.close
    @dst.unlink
  end
end
