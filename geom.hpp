#pragma once

namespace boom {
	//! calculates reciprocal value
	template <class T>
	struct Reciprocal {
		T operator()(const T& v) const {
			return 1 / v;
		}
	};
	template <>
	struct Reciprocal<float> {
		float operator()(float v) const {
			return spn::_sseRcp22Bit(v);
		}
	};
	//		template <int M, int N, bool Align>
	//		struct Reciprocal<MatT<M,N,Align>> {
	//			using MT = MatT<M,N,Align>;
	//			MT operator()(const MT& m) const {
	//				return m;
	//			}
	//		};

	//! stores arbitrary value and it's reciprocal value
	template <class T>
	struct RecipValue {
		T	value[2];

		RecipValue() = default;
		void operator = (const T& v) {
			value[0] = v;
			Reciprocal<T> recip;
			value[1] = recip(v);
		}
		const T& get(bool bInv) const {
			return value[static_cast<int>(bInv)];
		}
	};
	//! to gain compatibility with RecipValue
	template <class T>
	struct NmlValue {
		T	value;

		NmlValue() = default;
		NmlValue(const T& t): value(t) {}
		void operator = (const T& v) {
			value = v;
		}
		T get(bool bInv) const {
			if(!bInv)
				return value;
			Reciprocal<T> recip;
			return recip(value);
		}
		operator T& () {
			return value;
		}
	};

	template <class Core, template <class> class IM>
	class Cache : public IM<Core> {
		protected:
			using Info = typename Core::cache_type;
			Core				_core;
			mutable Info		_info;
			mutable uint32_t	_rflag;

		public:
			Cache(): _rflag(RFL_ALL) {}
			template <class... T>
			Cache(const T&... args): _core(args...), _rflag(RFL_ALL) {}
			Cache(const Cache& c) = default;
			Cache(const Core& src): _core(src), _rflag(RFL_ALL) {}

			const Core& getCore() const { return _core; }

			#define SEQ_CACHE ((getArea, GetArea, area, RFL_AREA, 1))\
							((getInertia, GetInertia, inertia, RFL_INERTIA, 1))\
							((getCenter, GetCenter, center, RFL_CENTER, 0))\
							((getBCircle, GetBCircle, bcircle, RFL_BCIRCLE, 0))
			#define DEF_CACHE2(func, stFunc, field, flag, arg)	decltype(_info.field.get(0)) func(BOOST_PP_IF(arg, bool bInv=false, NOTHING)) const override {}
//					if(_rflag & flag) spn::Bit::Clear(_rflag, flag | Core::stFunc(_info));
//					return _info.field[BOOST_PP_IF(arg, bInv, NOTHING)]; }
			#define DEF_CACHE(z,dummy,elem)	DEF_CACHE2(BOOST_PP_TUPLE_ELEM(0,elem),BOOST_PP_TUPLE_ELEM(1,elem),BOOST_PP_TUPLE_ELEM(2,elem),BOOST_PP_TUPLE_ELEM(3,elem),BOOST_PP_TUPLE_ELEM(4,elem))
			BOOST_PP_SEQ_FOR_EACH(DEF_CACHE, 0, SEQ_CACHE)
	};
}
