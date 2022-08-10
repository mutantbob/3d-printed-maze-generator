use crate::Point3D;

enum RingCommand {
    Append(Point3D),
    Unshift(Point3D),
    Close,
    Nope,
}

#[derive(Default)]
pub struct RingAccumulator {
    pub open_strings: Vec<Vec<Point3D>>,
    pub closed_strings: Vec<Vec<Point3D>>,
}

impl RingAccumulator {
    pub fn absorb(&mut self, p1: Point3D, p2: Point3D) {
        for idx in 0..self.open_strings.len() {
            let string = &mut self.open_strings[idx];
            let cmd = Self::command_for(&p1, &p2, string);
            match cmd {
                RingCommand::Append(p) => {
                    string.push(p);
                    self.consolidate_open_strings();
                    return;
                }
                RingCommand::Unshift(p) => {
                    string.insert(0, p);
                    self.consolidate_open_strings();
                    return;
                }
                RingCommand::Close => {
                    self.closed_strings.push(self.open_strings.remove(idx));
                    return;
                }
                RingCommand::Nope => {}
            }
        }

        self.open_strings.push(vec![p1, p2]);
    }

    fn consolidate_open_strings(&mut self) {
        for i in 0..self.open_strings.len() {
            let a = self.open_strings[i].first().unwrap();
            let b = self.open_strings[i].last().unwrap();
            for j in (i + 1)..self.open_strings.len() {
                let n = self.open_strings[j].first().unwrap();
                let o = self.open_strings[j].last().unwrap();

                if Self::close_enough(a, n) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.remove(0);
                    self.open_strings[i].splice(0..0, doomed.into_iter().rev());
                    self.consolidate_open_strings();
                    return;
                } else if Self::close_enough(b, n) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.remove(0);
                    self.open_strings[i].extend(doomed);
                    self.consolidate_open_strings();
                    return;
                } else if Self::close_enough(a, o) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.pop();
                    self.open_strings[i].splice(0..0, doomed);
                    self.consolidate_open_strings();
                    return;
                } else if Self::close_enough(b, o) {
                    let mut doomed = self.open_strings.remove(j);
                    doomed.pop();
                    self.open_strings[i].extend(doomed.into_iter().rev());
                    self.consolidate_open_strings();
                    return;
                }
            }
        }
    }

    fn command_for(p1: &Point3D, p2: &Point3D, string: &[Point3D]) -> RingCommand {
        let begin = string.first().unwrap();
        let end = string.last().unwrap();
        if Self::close_enough(begin, p1) {
            if Self::close_enough(end, p2) {
                RingCommand::Close
            } else {
                RingCommand::Unshift(*p2)
            }
        } else if Self::close_enough(end, p1) {
            if Self::close_enough(begin, p2) {
                RingCommand::Close
            } else {
                RingCommand::Append(*p2)
            }
        } else if Self::close_enough(begin, p2) {
            if Self::close_enough(end, p1) {
                RingCommand::Close
            } else {
                RingCommand::Unshift(*p1)
            }
        } else if Self::close_enough(end, p2) {
            if Self::close_enough(begin, p1) {
                RingCommand::Close
            } else {
                RingCommand::Append(*p1)
            }
        } else {
            RingCommand::Nope
        }
    }

    fn close_enough(p1: &Point3D, p2: &Point3D) -> bool {
        let delta = p1
            .iter()
            .zip(p2.iter())
            .map(|(a, b)| (a - b).abs())
            .fold(0.0, |a, b| if a < b { b } else { a });
        delta < 1.0e-6
    }
}
