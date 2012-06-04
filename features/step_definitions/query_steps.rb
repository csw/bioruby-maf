When /^filter for only the species$/ do |table|
  # table is a Cucumber::Ast::Table
  sp = table.raw.collect { |row| row[0] }
  @parser.sequence_filter = { :only_species => sp }
end
