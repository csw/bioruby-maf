source "http://rubygems.org"
# Add dependencies required to use your gem here.

gemspec

# Add dependencies to develop your gem here.
# Include everything needed to run rake, tests, features, etc.
group :development do
  gem "rdoc", "~> 3.12"
  gem "simplecov", "~> 0.6.4", :platforms => :mri
  gem "yard", "~> 0.8.1"
  gem "kramdown", "~> 0.13.6"
  gem "redcarpet", "~> 2.1.1", :platforms => :mri
  gem "ronn", "~> 0.7.3", :platforms => :mri
  gem "sinatra", "~> 1.3.2" # for ronn --server
end

group :test do
  gem "bundler", ">= 1.0.0"
  gem "rake", ">= 0.9"
  gem "cucumber", ">= 0"
  gem "rspec", "~> 2.10.0"
  gem "rubygems-tasks", "~> 0.2.3"
  gem "aruba", "~> 0.4.11"
end
