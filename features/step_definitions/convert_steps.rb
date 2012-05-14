Given /^a MAF source file "(.*?)"$/ do |src|
  @src = $test_data + src
  @src.exist?.should be_true
end

Given /^MAF data:$/ do |string|
  @src = string
end

When /^I select FASTA output$/ do
  @out_fmt = :fasta
end

When /^process the file$/ do
  @dst = Tempfile.new(['cuke', ".#{@out_fmt.to_s}"])
  proc = MAFProcessor.new
  # src can be a Pathname or a String with MAF data
  proc.src = @src
  proc.output_format = format
  proc.dst = @dst
  proc.process()
end

Then /^the output should match "(.*?)"$/ do |ref|
  ref_p = $test_data + ref
  ref.p.exist?.should be_true
  system("diff #{ref} #{@dst.path} >/dev/null 2>&1").should be_true
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
