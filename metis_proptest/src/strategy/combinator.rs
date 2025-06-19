use std::{marker::PhantomData, ptr::NonNull};

use crate::costs::CaseCost;

use super::{Case, Strategy, StrategyIter, StrategyKind};

pub struct DynStrategyIter<'a> {
    clone: unsafe fn (&DynStrategyIter<'a>) -> DynStrategyIter<'a>,
    data: NonNull<u8>,
    next: unsafe fn (&mut DynStrategyIter<'a>) -> Option<Case>,
    drop: unsafe fn (&mut DynStrategyIter),

    _phantom: PhantomData<&'a ()>,
}

impl<'a> DynStrategyIter<'a> {
    pub fn new<I: StrategyIter + Sized + 'a> (it: I) -> Self {
        let it = Box::new(it);
        unsafe fn clone<'a, I: Clone> (d: &DynStrategyIter<'a>) -> DynStrategyIter<'a> {
            let data = unsafe { (*(&raw const d.data).cast::<Box<I>>()).clone() };
            let data = unsafe { NonNull::new_unchecked(Box::into_raw(data).cast::<u8>()) };
            DynStrategyIter { clone: d.clone, data, next: d.next, drop: d.drop, _phantom: PhantomData }
        }
        unsafe fn next<'a, I: Iterator<Item = Case> + Sized> (d: &mut DynStrategyIter<'a>) -> Option<Case> {
            let data = unsafe { (d.data.cast::<I>()).as_mut() };
            data.next()
        }
        unsafe fn drop<'a, I: Sized> (d: &mut DynStrategyIter<'a>) {
            let data = (&raw mut d.data).cast::<Box<I>>() ;
            unsafe { std::ptr::drop_in_place(data) };
        }
        Self {
            clone: clone::<I>,
            data: unsafe { NonNull::new_unchecked(Box::into_raw(it)).cast() },
            next: next::<I>,
            drop: drop::<I>,
            _phantom: PhantomData,
        }
    }
}


impl<'a> Iterator for DynStrategyIter<'a> {
    type Item = Case;

    fn next(&mut self) -> Option<Self::Item> {
        unsafe { (self.next)(self) }
    }
}

impl<'a> Clone for DynStrategyIter<'a> {
    fn clone(&self) -> Self {
        unsafe { (self.clone)(self) }
    }
}

type OpaqueDynStrategyIter = [u8; size_of::<DynStrategyIter>()];

struct DynStrategyVTable {
    dont_backtrack: unsafe fn (NonNull<u8>) -> StrategyKind,
    reduction_power: unsafe fn (NonNull<u8>, case: &Case) -> CaseCost,
    is_valid: unsafe fn (NonNull<u8>, case: &Case) -> bool,
    apply: unsafe fn (NonNull<u8>, case: Case) -> Option<OpaqueDynStrategyIter>,
    drop_data: Option<unsafe fn(NonNull<u8>)>,
}

pub struct DynStrategy<'a> {
    _item: PhantomData<&'a ()>,
    data: NonNull<u8>,
    vtable: &'static DynStrategyVTable
}

impl<'a> DynStrategy<'a> {
    /// since `Box<S>` has the same ABI as `&S`, we just need to add drop and clone impls for boxed
    const fn construct_vtable_for_ref<S: Strategy + Sized>() -> DynStrategyVTable {
        unsafe fn cast_data<'a, S>(p: NonNull<u8>) -> &'a S {
            unsafe { p.cast::<S>().as_ref() }
        }

        unsafe fn dont_backtrack<S: Strategy + Sized>(p: NonNull<u8>) -> StrategyKind {
            unsafe { cast_data::<S>(p).kind() }
        }

        unsafe fn reduction_power<S: Strategy + Sized>(p: NonNull<u8>, case: &Case) -> CaseCost {
            unsafe { cast_data::<S>(p).cost(case) }
        }

        unsafe fn is_valid<S: Strategy + Sized>(p: NonNull<u8>, case: &Case) -> bool {
            unsafe { cast_data::<S>(p).is_valid(case) }
        }

        unsafe fn apply<'a, S: Strategy + Sized + 'a>(p: NonNull<u8>, case: Case) -> Option<OpaqueDynStrategyIter> {
            unsafe { cast_data::<S>(p).apply(case).map(DynStrategyIter::new).map(|it| std::mem::transmute(it)) }
        }

        DynStrategyVTable { 
            dont_backtrack: dont_backtrack::<S>,
            reduction_power: reduction_power::<S>,
            is_valid: is_valid::<S>,
            apply: apply::<S>,
            drop_data: None,
        }
    }

    pub const fn new<S: Strategy + Sized>(strat: &'a S) -> Self {
        let vt: &'static _ = & const { Self::construct_vtable_for_ref::<S>() };

        Self {
            _item: PhantomData,
            data: NonNull::from_ref(strat).cast(),
            vtable: vt,
        }
    }

    pub fn new_owned<S: Strategy + Sized + 'a>(strat: S) -> Self {
        Self::new_from_box(Box::new(strat))
    }

    pub fn new_from_box<S: Strategy + Sized + 'a>(strat: Box<S>) -> Self {
        unsafe fn drop_data<S>(p: NonNull<u8>) {
            let _ = unsafe { Box::from_raw(p.cast::<S>().as_ptr()) };
        }
        let vt: &'static _ = & const { 
            DynStrategyVTable {
                drop_data: Some(drop_data::<S>),
                ..Self::construct_vtable_for_ref::<S>()
            }
        };
        let data = unsafe { NonNull::new_unchecked(Box::into_raw(strat)) };
        Self {
            _item: PhantomData,
            data: data.cast(),
            vtable: vt,
        }
    }
}

impl Drop for DynStrategy<'_> {
    fn drop(&mut self) {
        if let Some(drop_fn) = self.vtable.drop_data {
            unsafe { drop_fn(self.data) }
        }
    }
}

impl<'a> Strategy for DynStrategy<'a> {
    fn cost(&self, case: &Case) -> CaseCost {
        unsafe { (self.vtable.reduction_power)(self.data, case) }
    }

    fn kind(&self) -> StrategyKind {
        unsafe { (self.vtable.dont_backtrack)(self.data) }
    }

    fn is_valid(&self, case: &Case) -> bool {
        unsafe { (self.vtable.is_valid)(self.data, case) }
    }

    #[allow(refining_impl_trait_internal)]
    fn apply(&self, case: Case) -> Option<DynStrategyIter<'a>> {
        unsafe { (self.vtable.apply)(self.data, case).map(|i| std::mem::transmute::<_, DynStrategyIter<'a>>(i)) }
    }
}

impl<'a, S: Strategy> From<&'a S> for DynStrategy<'a> {
    fn from(value: &'a S) -> Self {
        DynStrategy::new(value)
    }
}

impl<'a, S: Strategy + 'a> From<Box<S>> for DynStrategy<'a> {
    fn from(value: Box<S>) -> Self {
        DynStrategy::new_from_box(value)
    }
}

