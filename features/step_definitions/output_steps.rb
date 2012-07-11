When /^open a new MAF writer$/ do
  @dst = Tempfile.new(["cuke", ".maf"])
  @writer = Bio::MAF::Writer.new(@dst)
end

When /^write the header from the original MAF file$/ do
  @writer.write_header(@parser.header)
end

When /^write a default header$/ do
  @writer.write_header(Bio::MAF::Header.default)
end

When /^write all the parsed blocks$/ do
  @writer.write_blocks(@parser.parse_blocks)
end

When /^write all the matched blocks$/ do
  @writer.write_blocks(@blocks)
end

RSpec::Matchers.define :match_except_ws do |expected|
  match do |actual|
    system("diff --ignore-space-change --brief #{expected} #{actual} >/dev/null 2>&1")
  end

  failure_message_for_should do |actual|
    msg = "File contents did not match. Diff:\n"
    msg << `diff --unified --ignore-space-change #{expected} #{actual}`
  end
end

Then /^the output should match, except whitespace, "(.+)"$/ do |ref|
  @dst.path.should match_except_ws($test_data + ref)
end
