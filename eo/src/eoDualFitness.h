/*

(c) 2010 Thales group

    This library is free software; you can redistribute it and/or
    modify it under the terms of the GNU Lesser General Public
    License as published by the Free Software Foundation; version 2 
    of the License.

    This library is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    Lesser General Public License for more details.

    You should have received a copy of the GNU Lesser General Public
    License along with this library; if not, write to the Free Software
    Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

Contact: http://eodev.sourceforge.net

Authors:
    Johann Dréo <johann.dreo@thalesgroup.com>

*/

#ifndef _eoDualFitness_h_
#define _eoDualFitness_h_

#include <functional>
#include <iostream>
#include <utility> // for std::pair

//! A fitness class that permits to compare feasible and unfeasible individuals. It guaranties that a feasible individual will always be better than an unfeasible one.
/**
 * Use this class as fitness if you have some kind of individuals 
 * that must be always considered as better than others while having the same fitness type.
 *
 * Wraps a scalar fitness _values such as a double or int, with the option of
 * maximizing (using less<BaseType>, @see eoMaximizingDualFitness)
 * or minimizing (using greater<BaseType>, @see eoMinimizingDualFitness).
 *
 * Suitable constructors, assignments and casts are defined to work
 * with those quantities as if they were a pair of: a BaseType and a boolean.
 *
 * When changing the fitness, you can use:
 *     individual.fitness( std::make_pair<BaseType,bool>( fitness, feasibility ) );
 * 
 * Be aware that, when printing or reading an eDualFitness instance on a iostream,
 * friend IO classes use a space separator.
 *
 * This class overrides operator<() to use the Compare template argument and handle feasibility.
 * Over operators are coded using this sole function.
*/
template <class BaseType, class Compare >
class eoDualFitness
{
private:
    //! Scalar type of the fitness (generally a double)
    BaseType _value;

    //! Flag that marks if the individual is feasible
    bool _is_feasible;

public:

    //! Empty initialization
    /*!
     * Unfeasible by default
     */
    eoDualFitness() : 
        _value(), 
        _is_feasible(false) 
    {}

    //! Copy constructor
    eoDualFitness(const eoDualFitness& other) :
        _value(other._value),
        _is_feasible(other._is_feasible) 
    {}

    //! Constructor from explicit value/feasibility
    eoDualFitness(const BaseType& v, const bool& is_feasible) : 
        _value(v),
        _is_feasible(is_feasible)
    {}

    //! From a std::pair (first element is the value, second is the feasibility)
    eoDualFitness(const std::pair<BaseType,bool>& dual) :
        _value(dual.first),
        _is_feasible(dual.second)
    {}

    //! Copy operator
    eoDualFitness& operator=(const eoDualFitness& other)
    {
        _value = other._value;
        _is_feasible = other._is_feasible;
       return *this; 
    }

    //! Copy operator from a std::pair
    eoDualFitness& operator=(const std::pair<BaseType,bool>& v)
    {
        _value = v.first;
        _is_feasible = v.second;
        return *this; 
    }

    //! Comparison that separate feasible individuals from unfeasible ones. Feasible are always better
    /*!
     * Use less as a default comparison operator 
     * (see the "Compare" template of the class to change this behaviour,
     * @see eoMinimizingDualFitness for an example).
     */
    bool operator<(const eoDualFitness& other) const
    { 
        // am I better (less, by default) than the other ?

        // if I'm feasible and the other is not
        if( this->_is_feasible && !other._is_feasible ) {
            // no, the other has a better fitness
            return false;

        } else if( !this->_is_feasible && other._is_feasible ) {
            // yes, a feasible fitness is always better than an unfeasible one
            return true;

        } else { 
            // the two fitness are of the same type
            // lets rely on the comparator
            return Compare()(_value, other._value);
        }
    }

    //! Greater: if the other is lesser than me
    bool operator>( const eoDualFitness<BaseType, Compare>& other ) const  { return other < *this; }

    //! Less or equal: if the other is not lesser than me
    bool operator<=( const eoDualFitness<BaseType, Compare>& other ) const { return !(other < *this); }

    //! Greater or equal: if the other is not greater than me
    bool operator>=(const eoDualFitness<BaseType, Compare>& other ) const { return !(*this < other); }

public:

    //! Print an eoDualFitness instance as a pair of numbers, separated by a space
    template <class F, class Cmp>
    friend
    std::ostream& operator<<(std::ostream& os, const eoDualFitness<F, Cmp>& f)
    {
        os << f._value << " " << f._is_feasible;
        return os;
    }

    //! Read an eoDualFitness instance as a pair of numbers, separated by a space
    template <class F, class Cmp>
    friend 
    std::istream& operator>>(std::istream& is, eoDualFitness<F, Cmp>& f)
    {
        F value;
        is >> value;

        bool feasible;
        is >> feasible;

        f = std::make_pair<F,bool>( value, feasible );
        return is;
    }

};

//! Compare dual fitnesses as if we were maximizing
typedef eoDualFitness<double, std::less<double> >    eoMaximizingDualFitness;

//! Compare dual fitnesses as if we were minimizing
typedef eoDualFitness<double, std::greater<double> > eoMinimizingDualFitness;

#endif // _eoDualFitness_h_
