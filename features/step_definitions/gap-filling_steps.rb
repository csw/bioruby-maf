Given /^chromosome reference sequence:$/ do |string|
  @refseq = Bio::FastaFormat.new(string)
end

When /^tile ([^:\s]+):(\d+)-(\d+)( with the chromosome reference)?$/ do |seq, i_start, i_end, ref_p|
  @tiler = Bio::MAF::Tiler.new
  @tiler.index = @idx
  @tiler.parser = @parser
  @tiler.reference = @refseq.seq if ref_p
  @tiler.interval = Bio::GenomicInterval.zero_based(seq,
                                                    i_start.to_i,
                                                    i_end.to_i)
end

When /^tile with species \[(.+?)\]$/ do |species_text|
  @tiler.species = species_text.split(/,\s*/)
end

When /^map species (\S+) as (\S+)$/ do |sp1, sp2|
  @tiler.species_map[sp1] = sp2
end

When /^write the tiled data as FASTA$/ do
  @dst = Tempfile.new(["cuke", ".fa"])
  @tiler.write_fasta(@dst)
end

Then /^the FASTA data obtained should be:$/ do |string|
  @dst.seek(0)
  @dst.read.rstrip.should == string.rstrip
end
