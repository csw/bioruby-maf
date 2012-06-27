source "http://rubygems.org"
# Add dependencies required to use your gem here.
gem "bio-bigbio"
# temporarily using local copies
#gem "bio-ucsc-api", "~> 0.4.0"
gem "bio-genomic-interval", "~> 0.1.2"
gem "kyotocabinet-ruby", "~> 1.27.1", :platforms => [:mri, :rbx]
gem "kyotocabinet-java", "~> 0.2.0", :platforms => :jruby


# Add dependencies to develop your gem here.
# Include everything needed to run rake, tests, features, etc.
group :development do
  gem "rdoc", "~> 3.12"
  gem "simplecov", "~> 0.6.4", :platforms => :mri
  gem "yard", "~> 0.8.1"
  gem "kramdown", "~> 0.13.6"
  gem "ronn", "~> 0.7.3"
  gem "sinatra", "~> 1.3.2" # for ronn --server
end

group :test do
  gem "jeweler", "~> 1.8.3"
  gem "bundler", ">= 1.0.0"
  gem "rake", ">= 0.9"
  gem "cucumber", ">= 0"
  gem "rspec", "~> 2.10.0"
end
