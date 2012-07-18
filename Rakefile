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

require 'rubygems/tasks'
# we only want to do the SCM tag/push stuff once, on MRI
use_scm = (RUBY_PLATFORM != 'java')
Gem::Tasks.new(:scm => {:tag => use_scm, :push => use_scm})

require 'rspec/core'
require 'rspec/core/rake_task'
RSpec::Core::RakeTask.new(:spec) do |spec|
  spec.pattern = FileList['spec/**/*_spec.rb']
end

require 'cucumber/rake/task'
Cucumber::Rake::Task.new do |t|
  opts = "features"
  opts << ' --tags ~@no_jruby' if RUBY_PLATFORM == 'java'
  t.cucumber_opts = opts
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
    #sh "ronn --roff --html --style toc --date #{$gemspec.date.strftime('%Y-%m-%d')} --manual='BioRuby Manual' --organization='#{$gemspec.author}' #{file_spec}"
    sh "ronn --roff --html --style toc --date #{Time.now.strftime('%Y-%m-%d')} --manual='BioRuby Manual' --organization='BioRuby' #{file_spec}"
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
