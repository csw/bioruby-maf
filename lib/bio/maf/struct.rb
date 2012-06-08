## NOTE: this is probably not the best place for this, ultimately.
## If it works, think about moving it.

module Bio

  module MAF

    class Struct
      def initialize(spec)
        @members = []
        @by_name = {}
        offset = 0
        spec.each do |m_spec|
          m = Member.new(offset, *m_spec)
          @members << m
          @by_name[m.name] = m
          offset += m.size
        end
      end

      def extractor_fmt(*names)
        extract = names.collect { |name| @by_name.fetch(name) }
        extract.sort_by! { |m| m.offset }
        fmt = ''
        pos = 0
        extract.each do |member|
          if member.offset != pos
            fmt << "@#{member.offset}"
            pos = member.offset
          end
          fmt << member.fmt
          pos += member.size
        end
        return fmt
      end
    end

    TYPE_PROPS = {
      :uint8  => { :size => 1, :fmt => 'C'  },
      :uint16 => { :size => 2, :fmt => 'S>' },
      :uint32 => { :size => 4, :fmt => 'L>' },
      :uint64 => { :size => 8, :fmt => 'Q>' }
    }

    class Member
      attr_reader :offset, :name, :type, :size, :fmt
      def initialize(offset, name, type)
        @offset = offset
        @name = name
        @type = type
        props = TYPE_PROPS.fetch(type)
        @size = props.fetch(:size)
        @fmt = props.fetch(:fmt)
      end
    end

 end
  
end
