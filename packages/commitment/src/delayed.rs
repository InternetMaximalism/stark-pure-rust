use std::cell::RefCell;
use std::ops::Deref;
use std::rc::Rc;

#[derive(Clone)]
pub enum Thunk<'a, T: 'a> {
  NotYet(Rc<dyn Fn() -> T + 'a>),
  Memo(T),
}

#[derive(Clone)]
pub struct Delayed<'a, T: 'a> {
  thunk: RefCell<Box<Thunk<'a, T>>>,
}

impl<'a, T: 'a> Delayed<'a, T> {
  pub fn new<F>(f: F) -> Self
  where
    F: Fn() -> T + 'a,
  {
    Delayed {
      thunk: RefCell::new(Box::new(Thunk::NotYet(Rc::new(f)))),
    }
  }

  pub fn force(&self) {
    let thunk = &mut *self.thunk.borrow_mut();
    let val = match **thunk {
      Thunk::NotYet(ref invoke) => Box::new(Thunk::Memo(invoke())),
      Thunk::Memo(_) => return,
    };
    *thunk = val;
  }
}

impl<'a, T: 'a> Deref for Delayed<'a, T> {
  type Target = T;
  fn deref(&self) -> &T {
    self.force();
    let thunk = unsafe { self.thunk.as_ptr().as_ref().unwrap() };
    match **thunk {
      Thunk::Memo(ref v) => v,
      _ => unreachable!(),
    }
  }
}

#[macro_export]
macro_rules! lazily {
  ($($b:tt)+) => {
    Delayed::new(move || { $($b)+ })
  }
}
