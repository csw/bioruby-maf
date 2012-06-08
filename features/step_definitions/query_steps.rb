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

When /^filter for blocks with text size at (least|most) (\d+)$/ do |op, len|
  constraint = case op
               when 'least' then :min_size
               when 'most' then :max_size
               else raise "bad operator #{op}!"
               end
  @block_filter = { constraint => len.to_i}
end

When /^filter for blocks with text size between (\d+) and (\d+)$/ do |min, max|
  @block_filter = {
    :min_size => min.to_i,
    :max_size => max.to_i
  }
end
