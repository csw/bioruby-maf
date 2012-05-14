Given /^a MAF source file "(.*?)"$/ do |src|
  @src_p = $test_data + src
  @src_p.exist?.should be_true
end

When /^I select FASTA output$/ do
  @out_fmt = :fasta
end

When /^process the file$/ do
  @dst = Tempfile.new(['cuke', ".#{@out_fmt.to_s}"])
  @src_p.open do |src_f|
    proc = MAFProcessor.new
    proc.src = src_f
    proc.output_format = format
    proc.dst = @dst
    proc.process()
  end
end

Then /^the output should match "(.*?)"$/ do |ref|
  ref_p = $test_data + ref
  ref.p.exist?.should be_true
  system("diff #{ref} #{@dst.path} >/dev/null 2>&1").should be_true
end
