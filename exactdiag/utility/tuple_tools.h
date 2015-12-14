#pragma once

#include <tuple>
#include <typeinfo>
#include <type_traits>
#include <utility>
#include <algorithm>
#include <functional>

#define USE_BINOP_INSIDE


namespace aux
{

template<std::size_t _Idx, typename _Head, typename... _Tail>
struct _Elementwise_impl
    : public _Elementwise_impl<_Idx+1, _Tail...>
{
  const _Head & _head_impl;

  template <typename _TupleType>
      _Elementwise_impl(_TupleType&& t)
      : _Elementwise_impl<_Idx+1, _Tail...>(std::forward<_TupleType>(t))
      , _head_impl(std::get<_Idx>(t))
  {
  }

  //! ========== Unary ==========
  template <typename FuncTuple>
      auto apply_unary(FuncTuple&& func) const
      -> decltype(
          std::tuple_cat(
            std::make_tuple(std::get<_Idx>(func)(_head_impl)),
            _Elementwise_impl<_Idx+1, _Tail...>::apply_unary(std::forward<FuncTuple>(func))
            )
        )
  {
    return std::tuple_cat(
        std::make_tuple(std::get<_Idx>(func)(_head_impl)),
        _Elementwise_impl<_Idx+1, _Tail...>::apply_unary(std::forward<FuncTuple>(func))
      );
  }


  //! ========== Binary ==========
  template <typename FuncTuple>
      auto apply_binary(FuncTuple && func, _Elementwise_impl const & rhs) const
      -> decltype(
          std::tuple_cat(
              std::make_tuple(std::get<_Idx>(func)(_head_impl, rhs._head_impl)),
              _Elementwise_impl<_Idx+1, _Tail...>::apply_binary(
                std::forward<FuncTuple>(func), 
                static_cast<_Elementwise_impl<_Idx + 1, _Tail...> const &>(rhs))
            )
        )
  {
    return std::tuple_cat(
        std::make_tuple(std::get<_Idx>(func)(_head_impl, rhs._head_impl)),
        _Elementwise_impl<_Idx+1, _Tail...>::apply_binary(
          std::forward<FuncTuple>(func),
          static_cast<_Elementwise_impl<_Idx + 1, _Tail...> const &>(rhs))
      );
  }

  bool all() const {
    return bool(_head_impl) && _Elementwise_impl<_Idx+1, _Tail...>::all();
  }

  bool any() const {
    return bool(_head_impl) || _Elementwise_impl<_Idx+1, _Tail...>::any();
  }
};

template<std::size_t _Idx, typename _Head>
struct _Elementwise_impl<_Idx, _Head>
{
  const _Head & _head_impl;

  template <typename _TupleType>
      _Elementwise_impl(_TupleType&& t)
      : _head_impl(std::get<_Idx>(t))
  {  }

  template <typename FuncTuple>
  auto apply_unary(FuncTuple&& func) const
      -> decltype(std::make_tuple(std::get<_Idx>(func)(_head_impl)))
  { return std::make_tuple(std::get<_Idx>(func)(_head_impl)); }

  template <typename FuncTuple>
  auto apply_binary(FuncTuple && func, _Elementwise_impl const & rhs) const
      -> decltype(std::make_tuple(std::get<_Idx>(func)(_head_impl, rhs._head_impl)))
  { return std::make_tuple(std::get<_Idx>(func)(_head_impl, rhs._head_impl)); }

  bool all() const { return bool(_head_impl); }
  bool any() const { return bool(_head_impl); }

};

} // namespace aux

template <typename ...Types>
using Elementwise = aux::_Elementwise_impl<0, Types...>;

//!
//! Binary Operations
//!
#define DEF_BINOP(OP, FUNC)                                             \
  template<typename ...Types>                                           \
  auto operator OP(const Elementwise<Types...> & lhs,                   \
                   const Elementwise<Types...> & rhs)                   \
      -> decltype( lhs.apply_binary( std::make_tuple(FUNC<Types>()...), rhs) ) \
  { return lhs.apply_binary( std::make_tuple(FUNC<Types>()...), rhs); }

DEF_BINOP(+ , std::plus)
DEF_BINOP(- , std::minus)
DEF_BINOP(* , std::multiplies)
DEF_BINOP(/ , std::divides)
DEF_BINOP(% , std::modulus)
DEF_BINOP(==, std::equal_to)
DEF_BINOP(!=, std::not_equal_to)
DEF_BINOP(<=, std::less_equal)
DEF_BINOP(>=, std::greater_equal)
DEF_BINOP(< , std::less)
DEF_BINOP(> , std::greater)
DEF_BINOP(&&, std::logical_and)
DEF_BINOP(||, std::logical_or)

#undef DEF_BINOP


template<typename ...Types>
auto elmax(const Elementwise<Types...> & lhs,
           const Elementwise<Types...> & rhs)
    -> decltype(lhs.apply_binary(
        std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::max)...),
        rhs))
{
  return lhs.apply_binary(
          std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::max)...),
          rhs);
}

template<typename ...Types>
auto elmin(const Elementwise<Types...> & lhs,
           const Elementwise<Types...> & rhs)
    -> decltype( lhs.apply_binary(
        std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::min)...),
        rhs) )
{
  return lhs.apply_binary(
          std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::min)...),
          rhs);
}

template<typename ...Types>
auto elmax(const std::tuple<Types...> & lhs,
           const std::tuple<Types...> & rhs)
    -> decltype( elementwise(lhs).apply_binary(
        std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::max)...),
        elementwise(rhs)) )
{
  return elementwise(lhs).apply_binary(
          std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::max)...),
          elementwise(rhs));
}


template<typename ...Types>
auto elmin(const std::tuple<Types...> & lhs,
           const std::tuple<Types...> & rhs)
    -> decltype( elementwise(lhs).apply_binary(
        std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::min)...),
        elementwise(rhs)) )
{
  return elementwise(lhs).apply_binary(
          std::make_tuple(static_cast<const Types&(*)(const Types&, const Types&) >(std::min)...),
          elementwise(rhs));
}



//!
//! Unary Operations
//!
#define DEF_UNIOP(OP, FUNC)                                             \
  template<typename ...Types>                                           \
  auto operator OP(const Elementwise<Types...> & lhs)                   \
      -> decltype( lhs.apply_unary( std::make_tuple(FUNC<Types>()...)) ) \
  { return lhs.apply_unary( std::make_tuple(FUNC<Types>()...)); }

DEF_UNIOP(-, std::negate)
DEF_UNIOP(!, std::logical_not)

#undef DEF_UNIOP


template <typename T, typename ...TR>
Elementwise<T, TR...>
elementwise(const std::tuple<T, TR...>& t)
{
  return Elementwise<T, TR...>(t);
}

template <typename ...Types>
bool all(const std::tuple<Types...>& t)
{
  return elementwise(t).all();
}
template <typename ...Types>
bool any(const std::tuple<Types...>& t)
{
  return elementwise(t).any();
}

