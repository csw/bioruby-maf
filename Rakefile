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
Jeweler::Tasks.new do |gem|
  # gem is a Gem::Specification... see http://docs.rubygems.org/read/chapter/20 for more options
  gem.name = "bio-maf"
  gem.homepage = "http://github.com/csw/bioruby-maf"
  gem.license = "MIT"
  gem.summary = %Q{TODO: one-line summary of your gem}
  gem.description = %Q{TODO: longer description of your gem}
  gem.email = "cswh@umich.edu"
  gem.authors = ["csw"]
  # dependencies defined in Gemfile
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

GEMSPEC = Gem::Specification.load('bio-maf.gemspec')

def ronn(format, file)
  sh "ronn --build #{format} --style toc --date #{GEMSPEC.date.strftime('%Y-%m-%d')} --manual='RubyGems Manual' --organization='#{GEMSPEC.author}' #{file}"
end

def ronn_files
  GEMSPEC.files.grep /^man\/.*\.ronn$/
end

def man_files
  ronn_files.map { |path| path.sub(/\.ronn$/, '') }
end

def html_man_files
  ronn_files.map { |path| path.sub(/\.ronn$/, '.html') }
end

rule %r{\.\d$} => "%p.ronn" do |task|
  ronn('--roff', task.source)
end

rule %r{\.\d.html$} => "%X.ronn" do |task|
  ronn('--html', task.source)
end

task :man => (man_files + html_man_files)

namespace :man do
  task :publish do
    html_man_files.each do |man|
      cp man, "../octopress/source/man/#{File.basename(man)}"
    end
  end
end
task 'man:publish' => :man

namespace :ronn do
  task :server do
    sh "ronn --server #{ronn_files.join(' ')}"
  end
end

#### RDoc (not currently used)

require 'rdoc/task'
Rake::RDocTask.new do |rdoc|
  version = File.exist?('VERSION') ? File.read('VERSION') : ""

  rdoc.rdoc_dir = 'rdoc'
  rdoc.title = "bio-maf #{version}"
  rdoc.rdoc_files.include('README*')
  rdoc.rdoc_files.include('lib/**/*.rb')
end
