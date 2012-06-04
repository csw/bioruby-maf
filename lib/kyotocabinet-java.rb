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

begin
  Java::Kyotocabinet::DB
rescue NameError
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

module Java::Kyotocabinet
  class Cursor
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
  end # class Cursor

  class DB
    alias_method :_match_prefix, :match_prefix
    def match_prefix(prefix, limit=-1)
      _match_prefix(prefix, limit)
    end

    def [](key)
      get(key)
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
