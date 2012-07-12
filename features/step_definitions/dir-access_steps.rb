Given /^indexed MAF files in "(.*?)"$/ do |dir|
  @opts ||= {}
  @access = Bio::MAF::Access.maf_dir(dir, @opts)
end

When /^I query for the genomic intervals$/ do |table|
  # table is a Cucumber::Ast::Table
  intervals = table.hashes.collect do |row|
    Bio::GenomicInterval.zero_based(row['chrom'],
                                    row['start'].to_i,
                                    row['end'].to_i)
  end
  @access.block_filter = @block_filter
  @blocks = @access.find(intervals).to_a
end
