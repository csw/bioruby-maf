##
## JRuby adaptation layer for Kyoto Cabinet Java interface.
##
## This requires the kyotocabinet-java library to be installed. This
## library builds a jar file, which should be on the Ruby $LOAD_PATH,
## the Java CLASSPATH, or in /usr/local/lib/ where it is placed by the
## default 'make install' target in kyotocabinet-java.
##
## kyotocabinet-java also builds a necessary JNI library,
## libjkyotocabinet.{jnilib,dylib}. This cannot be added to
## java.library.path at runtime, so either its installation path must
## be added to java.library.path with a JVM option like
## -Djava.library.path=/some/dir, or it must be copied to somewhere
## already on java.library.path.
##

unless RUBY_PLATFORM == 'java'
  raise "This library requires JRuby!"
end

require 'java'

unless Java.const_defined? :Kyotocabinet
  # need to require kyotocabinet.jar
  if File.exist? '/usr/local/lib/kyotocabinet.jar'
    require '/usr/local/lib/kyotocabinet.jar'
    begin
      Java::Kyotocabinet::DB
    rescue NameError
      raise "Kyoto Cabinet classes could not be loaded, probably due to JNI library linking problems. Verify that libjkyotocabinet is in a directory on java.library.path (#{java.lang.System.getProperty('java.library.path')})."
    end
  else
    raise "kyotocabinet.jar could not be found!"
  end
end

module KyotoCabinet
  module Adaptation
    BYTE_ARRAY = [1].to_java(:byte).java_class

    def self.included(mod)
      super(mod)
      mod.class_eval do
        def self.with_method_handle(mname, *args)
          $stderr.puts "in #{self}.with_method_handle(#{mname.inspect}, #{args.inspect})"
          mh = self.java_method(mname, args)
          $stderr.puts "got method handle #{mh}"
          yield mh
        end
      end
    end

  end
end

module Java::Kyotocabinet
  class Cursor
    include KyotoCabinet::Adaptation

    alias_method :_get, :get
    def get(step=false)
      r = self._get(step)
      if r
        return [String.from_java_bytes(r[0]),
                String.from_java_bytes(r[1])]
      else
        return nil
      end
    end

    alias_method :_get_key, :get_key
    def get_key(step=false)
      r = self._get_key(step)
      if r
        return String.from_java_bytes(r)
      else
        return nil
      end
    end

    alias_method :_jump, :jump
    def jump(key)
      self._jump(key.to_java_bytes)
    end
  end # class Cursor

  class DB
    include KyotoCabinet::Adaptation

    alias_method :_match_prefix, :match_prefix
    def match_prefix(prefix, limit=-1)
      self._match_prefix(prefix, limit)
    end

    alias_method :_get, :get
    def get(key)
      String.from_java_bytes(self._get(key.to_java_bytes))
    end
    alias_method :[], :get

    alias_method :_set, :set
    def set(k, v)
      self._set(k.to_java_bytes, v.to_s.to_java_bytes)
    end
    alias_method :[]=, :set

    alias_method :_set_bulk, :set_bulk
    def set_bulk(rec_h, atomic)
      ba = Java::byte[rec_h.size * 2, 0].new
      i = 0
      rec_h.each_pair do |k, v|
        ba[i]   = k.to_java_bytes
        ba[i+1] = v.to_java_bytes
        i += 2
      end
      self._set_bulk(ba, atomic)
    end

    alias_method :_synchronize, :synchronize
    def synchronize(hard=false, proc=nil)
      self._synchronize(hard, proc)
    end

    def cursor_process
      cur = self.cursor()
      begin
        yield cur
      ensure
        cur.disable()
      end
    end
  end # class DB
end

module KyotoCabinet
  import Java::Kyotocabinet
end
