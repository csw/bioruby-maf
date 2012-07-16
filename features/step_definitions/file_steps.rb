Given /^test files:$/ do |table|
  Pathname.new("tmp/aruba").mkpath
  table.raw.collect { |row| $test_data + row[0] }.each do |path|
    $stderr.puts "staging #{path}"
    system("cp #{path} tmp/aruba/")
  end
end
