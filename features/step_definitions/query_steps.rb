When /^filter for only the species$/ do |table|
  # table is a Cucumber::Ast::Table
  sp = table.raw.collect { |row| row[0] }
  @parser.sequence_filter = { :only_species => sp }
end

When /^filter for blocks with the species$/ do |table|
  # table is a Cucumber::Ast::Table
  sp = table.raw.collect { |row| row[0] }
  @block_filter = { :with_all_species => sp }
end

When /^filter for blocks with at least (\d+) sequences$/ do |n|
  @block_filter = { :at_least_n_sequences => n.to_i }
end
