

#ifndef VOTCA_CSG_TOPOLOGYOBJECTFACTORY
#define VOTCA_CSG_TOPOLOGYOBJECTFACTORY

#include <iostream>
#include <map>
#include <string>
#include <vector>

#include <boost/lexical_cast.hpp>

namespace votca {

namespace csg {
template <class key_t, class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
class ObjectFactory {
 private:
  typedef T<Bead_T, Molecule_T, Topology_T> *(*creator_t)();

 public:
  typedef T<Bead_T, Molecule_T, Topology_T> abstract_type;
  ObjectFactory(){};
  ~ObjectFactory(){};

  T<Bead_T, Molecule_T, Topology_T> *Create(const key_t &key);

  bool IsRegistered(const key_t &id) const;

  void Register(const key_t &key, creator_t creator);

  // template< class obj_t >
  template <template <class, template <class> class,
                      template <class, template <class> class> class>
            class obj_t>
  void Register(const key_t &key);

  static ObjectFactory<key_t, Bead_T, Molecule_T, Topology_T, T> &Instance() {
    static ObjectFactory<key_t, Bead_T, Molecule_T, Topology_T, T> this_;
    return this_;
  }

  const std::map<key_t, creator_t> &getObjects() { return objects_; }

 private:
  std::map<key_t, creator_t> objects_;
};

template <class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T,
          template <class, template <class> class,
                    template <class, template <class> class> class>
          class parent>
parent<Bead_T, Molecule_T, Topology_T> *create_policy_new() {
  return *new T<Bead_T, Molecule_T, Topology_T>();
}

template <class key_t, class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
inline T<Bead_T, Molecule_T, Topology_T>
    *ObjectFactory<key_t, Bead_T, Molecule_T, Topology_T, T>::Create(
        const key_t &key) {
  typename std::map<key_t, creator_t>::const_iterator it(objects_.find(key));
  if (it != objects_.end())
    return (it->second)();
  else
    throw std::runtime_error(
        "factory key " + boost::lexical_cast<std::string>(key) + " not found.");
}

template <class key_t, class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
inline bool ObjectFactory<key_t, Bead_T, Molecule_T, Topology_T,
                          T>::IsRegistered(const key_t &id_) const {
  return (objects_.find(id_) != objects_.end());
}

template <class key_t, class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
inline void ObjectFactory<key_t, Bead_T, Molecule_T, Topology_T, T>::Register(
    const key_t &key, creator_t creator) {
  (void)objects_
      .insert(typename std::map<key_t, creator_t>::value_type(key, creator))
      .second;
}

template <class key_t, class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
template <template <class, template <class> class,
                    template <class, template <class> class> class> class obj_t>
void ObjectFactory<key_t, Bead_T, Molecule_T, Topology_T, T>::Register(
    const key_t &key) {
  Register(key, create_policy_new<obj_t<Bead_T, Molecule_T, Topology_T>>);
}

template <class Bead_T, template <class> class Molecule_T,
          template <class, template <class> class> class Topology_T,
          template <class, template <class> class,
                    template <class, template <class> class> class> class T>
class ObjectFactoryRegister {
 public:
  template <class factory_type, class key_type>
  ObjectFactoryRegister(factory_type &factory, key_type &key) {
    factory.Register(key,
                     &create_policy_new<typename factory_type::abstract_type,
                                        T<Bead_T, Molecule_T, Topology_T>>);
  }
};

}  // namespace csg
}  // namespace votca

#endif  // VOTCA_CSG_TOPOLOGYOBJECTFACTORY
