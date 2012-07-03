Given /^chromosome reference sequence:$/ do |string|
  @refseq = Bio::FastaFormat.new(string)
end

When /^tile ([^:\s]+):(\d+)-(\d+)( with the chromosome reference)?$/ do |seq, i_start, i_end, ref_p|
  @tiler = Bio::MAF::Tiler.new
  @tiler.index = @index
  @tiler.parser = @parser
  @tiler.reference = @refseq if ref_p
  @tiler.interval = Bio::GenomicInterval.zero_based(seq,
                                               i_start.to_i,
                                               i_end.to_i)
end

When /^write the tiled data as FASTA$/ do
  pending
  @dst = Tempfile.new(["cuke", ".fa"])
  @tiler.write_fasta(@dst)
end

Then /^the FASTA data obtained should be:$/ do |string|
  pending
end
