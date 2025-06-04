#     Copyright (C) 2023 Jordi Llopis-Lorente. Contact: jorllolo@etsii.upv.es
#
#     This program is distributed in the hope that it will be useful,
#     but WITHOUT ANY WARRANTY; without even the implied warranty of
#     MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#     GNU General Public License for more details.
#
#     You should have received a copy of the GNU General Public License
#     along with this program.  If not, see <https://www.gnu.org/licenses/>

getQnet <- function(time, IsJs) {
  # qnet_mean: q_net value across the beats saved
  # qnet_beats: vector containing the q_net value for each of the beats saved
  # time: time vector
  # IsJs: matrix containing the currents during the beats saved

  qnet_beats <- numeric(length(IsJs))

  for (i in seq_along(IsJs)) {
    t <- time[[i]]
    Inet <- rowSums(IsJs[[i]][, c(2, 3, 4, 7, 8, 9)])
    qnet_total <- cumsum(diff(t) / 1000 * head(Inet, -1)) # NB: cumtrapz does an automatic trapezoidal rule, while cumsum() in R is a basic summation.
    qnet_beats[i] <- tail(qnet_total, 1)
  }

  qnet_mean <- mean(qnet_beats, na.rm = TRUE)

  return(list(qnet_mean = qnet_mean, qnet_beats = qnet_beats))
}
