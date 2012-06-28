# encoding: utf-8

require 'rubygems'
require 'bundler'
begin
  Bundler.setup(:default, :development)
rescue Bundler::BundlerError => e
  $stderr.puts e.message
  $stderr.puts "Run `bundle install` to install missing gems"
  exit e.status_code
end
require 'rake'

require 'jeweler'
$gemspec = nil
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "bio-maf"
  gem.homepage = "http://github.com/csw/bioruby-maf"
  gem.license = "MIT"
  gem.summary = %Q{MAF parser for BioRuby}
  gem.description = %Q{Multiple Alignment Format parser for BioRuby.}
  gem.email = "cswh@umich.edu"
  gem.authors = ["Clayton Wheeler"]
  # dependencies defined in Gemfile
  # kludgy, but it works
  $gemspec = gem
end
Jeweler::RubygemsDotOrgTasks.new

require 'rspec/core'
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end

require 'cucumber/rake/task'
Cucumber::Rake::Task.new do |features|
end

task :test => [ :spec, :cucumber ] 
task :default => :test

#### Man pages
# (borrowed from matthewtodd/shoe)
ronn_avail = begin
               require 'ronn'
               true
             rescue LoadError
               false
             end

if ronn_avail
  RONN_FILES = Rake::FileList["man/*.?.ronn"]

  desc "Generate man pages"
  task :man do
    file_spec = RONN_FILES.join(' ')
    sh "ronn --roff --html --style toc --date #{$gemspec.date.strftime('%Y-%m-%d')} --manual='BioRuby Manual' --organization='#{$gemspec.author}' #{file_spec}"
  end

  namespace :man do
    desc "Publish man pages to Octopress source dir"
    task :publish do
      RONN_FILES.map { |path| path.sub(/\.ronn$/, '.html') }.each do |man|
        cp man, "../octopress/source/man/#{File.basename(man)}"
      end
    end
  end
  task 'man:publish' => :man

  namespace :ronn do
    task :server do
      sh "ronn --server #{RONN_FILES.join(' ')}"
    end
  end
end # if ronn_avail

#### RDoc (not currently used)

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "bio-maf #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
