test_that("load_vmdb", {
  testthat::skip_if_not_installed("RSQLite")
  con <- DBI::dbConnect(RSQLite::SQLite(), ":memory:")
  on.exit(DBI::dbDisconnect(con), add = TRUE)

  # Minimal schema compatible with queries used in load_vmdb_*
  DBI::dbExecute(con, '
    CREATE TABLE obs_session (
        id integer PRIMARY KEY,
        longitude real NOT NULL,
        latitude real NOT NULL,
        elevation real NOT NULL,
        observer_id integer NULL,
        observer_name TEXT NULL,
        country TEXT NOT NULL,
        city TEXT NOT NULL
  )')
  DBI::dbExecute(con, '
    CREATE TABLE rate (
      id integer NOT NULL,
      shower varchar(6) NULL,
      period_start timestamp NOT NULL,
      period_end timestamp NOT NULL,
      sl_start double precision NOT NULL,
      sl_end double precision NOT NULL,
      session_id integer NOT NULL,
      freq integer NOT NULL,
      lim_mag real NOT NULL,
      t_eff real NOT NULL,
      f real NOT NULL,
      sidereal_time double precision NOT NULL,
      sun_alt double precision NOT NULL,
      sun_az double precision NOT NULL,
      moon_alt double precision NOT NULL,
      moon_az double precision NOT NULL,
      moon_illum double precision NOT NULL,
      field_alt double precision NULL,
      field_az double precision NULL,
      rad_alt double precision NULL,
      rad_az double precision NULL,
      CONSTRAINT rate_pkey PRIMARY KEY (id),
      CONSTRAINT rate_session_fk FOREIGN KEY (session_id)
          REFERENCES obs_session(id) MATCH SIMPLE
          ON UPDATE CASCADE
          ON DELETE CASCADE
  )')
  DBI::dbExecute(con, '
    CREATE TABLE magnitude (
      id integer NOT NULL,
      shower varchar(6) NULL,
      period_start timestamp NOT NULL,
      period_end timestamp NOT NULL,
      sl_start double precision NOT NULL,
      sl_end double precision NOT NULL,
      session_id integer NOT NULL,
      freq integer NOT NULL,
      mean double precision NOT NULL,
      lim_mag real NULL,
      CONSTRAINT magnitude_pkey PRIMARY KEY (id),
      CONSTRAINT magnitude_session_fk FOREIGN KEY (session_id)
          REFERENCES obs_session(id) MATCH SIMPLE
          ON UPDATE CASCADE
          ON DELETE CASCADE
  )')
  DBI::dbExecute(con, '
    CREATE TABLE magnitude_detail (
    id integer NOT NULL,
    magn integer NOT NULL,
    freq real NOT NULL,
    CONSTRAINT magnitude_detail_pkey PRIMARY KEY (id, magn),
    CONSTRAINT magnitude_detail_fk FOREIGN KEY (id)
        REFERENCES magnitude(id) MATCH SIMPLE
        ON UPDATE CASCADE
        ON DELETE CASCADE
  )')
  DBI::dbExecute(con, '
    CREATE TABLE rate_magnitude (
      rate_id integer NOT NULL,
      magn_id integer NOT NULL,
      "equals" boolean NOT NULL,
      CONSTRAINT rate_magnitude_pkey PRIMARY KEY (rate_id),
      CONSTRAINT rate_magnitude_rate_fk FOREIGN KEY (rate_id)
          REFERENCES rate (id) MATCH SIMPLE
          ON UPDATE CASCADE
          ON DELETE CASCADE,
      CONSTRAINT rate_magnitude_magn_fk FOREIGN KEY (magn_id)
          REFERENCES magnitude(id) MATCH SIMPLE
          ON UPDATE CASCADE
          ON DELETE CASCADE
  )')

  # Seed data: multiple sessions and observations to exercise filters
  DBI::dbExecute(con, "INSERT INTO obs_session (id, longitude, latitude, elevation, country, city, observer_id, observer_name)
    VALUES (1, 10.0, 50.0, 0.3, 'DE', 'Somewhere', 'XX', 'Doe, J.')")
  DBI::dbExecute(con, "INSERT INTO obs_session (id, longitude, latitude, elevation, country, city, observer_id, observer_name)
    VALUES (2, 11.0, 51.0, 0.5, 'DE', 'Elsewhere', 'YY', 'Roe, A.')")
  DBI::dbExecute(con, "INSERT INTO obs_session (id, longitude, latitude, elevation, country, city, observer_id, observer_name)
    VALUES (3, 12.0, 52.0, 0.1, 'FR', 'Paris', 'ZZ', 'Mete, O.')")
  DBI::dbExecute(con, "INSERT INTO rate (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, lim_mag,
      t_eff, f, sidereal_time, sun_alt, sun_az, moon_alt, moon_az, moon_illum, field_alt, field_az, rad_alt, rad_az)
    VALUES (100, 'PER', '2015-08-12 00:00:00', '2015-08-12 01:00:00', 140.0, 140.5, 1, 12, 6.5,
      0.9, 1.0, 10.0, -18.0, 120.0, -30.0, 200.0, 0.5, 50.0, 180.0, 45.0, 90.0)")
  DBI::dbExecute(con, "INSERT INTO rate (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, lim_mag,
      t_eff, f, sidereal_time, sun_alt, sun_az, moon_alt, moon_az, moon_illum, field_alt, field_az, rad_alt, rad_az)
    VALUES (101, 'PER', '2015-08-12 01:00:00', '2015-08-12 02:00:00', 140.6, 141.0, 2, 14, 5.8,
      0.8, 0.9, 11.0, -12.0, 122.0, -15.0, 201.0, 0.4, 48.0, 182.0, 44.5, 91.0)")
  DBI::dbExecute(con, "INSERT INTO rate (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, lim_mag,
      t_eff, f, sidereal_time, sun_alt, sun_az, moon_alt, moon_az, moon_illum, field_alt, field_az, rad_alt, rad_az)
    VALUES (102, 'GEM', '2015-12-14 00:00:00', '2015-12-14 01:00:00', 250.0, 250.4, 2, 20, 6.8,
      0.85, 1.1, 30.0, -6.0, 60.0, 12.0, 210.0, 0.8, 60.0, 175.0, 55.0, 88.0)")
  DBI::dbExecute(con, "INSERT INTO rate (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, lim_mag,
      t_eff, f, sidereal_time, sun_alt, sun_az, moon_alt, moon_az, moon_illum, field_alt, field_az, rad_alt, rad_az)
    VALUES (103, NULL, '2015-12-14 01:00:00', '2015-12-14 02:00:00', 350.0, 5.0, 3, 6, 6.1,
      0.7, 0.95, 31.0, -2.0, 61.0, -5.0, 211.0, 0.9, 58.0, 176.0, 15.0, 10.0)")
  DBI::dbExecute(con, 'INSERT INTO rate_magnitude (rate_id, magn_id, "equals") VALUES (100, 200, true)')
  DBI::dbExecute(con, 'INSERT INTO rate_magnitude (rate_id, magn_id, "equals") VALUES (101, 201, true)')
  DBI::dbExecute(con, 'INSERT INTO rate_magnitude (rate_id, magn_id, "equals") VALUES (102, 202, true)')
  DBI::dbExecute(con, "INSERT INTO magnitude (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, mean, lim_mag)
    VALUES (200, 'PER', '2015-08-12 00:00:00', '2015-08-12 01:00:00', 140.0, 140.5, 1, 12, 2.5, 6.5)")
  DBI::dbExecute(con, "INSERT INTO magnitude (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, mean, lim_mag)
    VALUES (201, 'PER', '2015-08-12 01:00:00', '2015-08-12 02:00:00', 140.6, 141.0, 2, 14, 2.7, 5.8)")
  DBI::dbExecute(con, "INSERT INTO magnitude (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, mean, lim_mag)
    VALUES (202, 'GEM', '2015-12-14 00:00:00', '2015-12-14 01:00:00', 250.0, 250.4, 2, 20, 3.2, 6.8)")
  DBI::dbExecute(con, "INSERT INTO magnitude (id, shower, period_start, period_end, sl_start, sl_end, session_id, freq, mean, lim_mag)
    VALUES (203, NULL, '2015-12-14 01:00:00', '2015-12-14 02:00:00', 350.0, 5.0, 3, 6, 2.0, 6.1)")
  DBI::dbExecute(con, "INSERT INTO magnitude_detail (id, magn, freq) VALUES
    (200, 3, 2.5), (200, 2, 4.0), (200, 1, 5.5),
    (201, 3, 3.0), (201, 2, 5.0), (201, 1, 6.0),
    (202, 4, 4.0), (202, 3, 8.0), (202, 2, 8.0),
    (203, 2, 4.0), (203, 1, 2.0)")

  # load_vmdb_rates: without extras
  res <- load_vmdb_rates(con, shower = 'PER')
  expect_type(res, 'list')
  expect_true(is.data.frame(res$observations))
  expect_null(res$sessions)
  expect_null(res$magnitudes)
  expect_equal(nrow(res$observations), 2)
  expect_setequal(res$observations$rate.id, c(100, 101))

  # with sessions and magnitudes
  res <- load_vmdb_rates(con, shower = 'PER', withSessions = TRUE, withMagnitudes = TRUE)
  expect_true(is.data.frame(res$sessions))
  expect_true(methods::is(res$magnitudes, 'table'))
  expect_equal(sort(row.names(res$sessions)), c('1', '2'))
  expect_equal(sort(row.names(res$magnitudes)), c('200', '201'))
  expect_setequal(colnames(res$magnitudes), c('1', '2', '3'))

  # shower filters (including sporadic via NA)
  expect_setequal(
    load_vmdb_rates(con, shower = c('PER', 'GEM'))$observations$rate.id,
    c(100, 101, 102)
  )
  expect_equal(
    load_vmdb_rates(con, shower = c(NA))$observations$rate.id,
    c(103)
  )
  expect_setequal(
    load_vmdb_rates(con, shower = c('PER', NA))$observations$rate.id,
    c(100, 101, 103)
  )

  # period filters
  expect_equal(
    load_vmdb_rates(con, period = c('2015-08-12 00:00:00', '2015-08-12 01:10:00'))$observations$rate.id,
    c(100)
  )
  expect_setequal(
    load_vmdb_rates(
      con,
      period = c(
        '2015-08-12 00:00:00', '2015-08-12 02:00:00',
        '2015-12-14 00:00:00', '2015-12-14 02:00:00'
      )
    )$observations$rate.id,
    c(100, 101, 102, 103)
  )

  # sl filters
  expect_equal(
    load_vmdb_rates(con, shower = 'PER', sl = c(140.0, 141.1))$observations$rate.id,
    c(101)
  )
  expect_equal(
    load_vmdb_rates(con, shower = c(NA), sl = c(340.0, 10.0))$observations$rate.id,
    c(103)
  )

  # limiting magnitude filters
  expect_setequal(
    load_vmdb_rates(con, lim.magn = c(6.0, 6.6))$observations$rate.id,
    c(100, 103)
  )

  # solar and lunar altitude filters
  expect_equal(
    load_vmdb_rates(con, sun.alt.max = -16)$observations$rate.id,
    c(100)
  )
  expect_setequal(
    load_vmdb_rates(con, moon.alt.max = 0)$observations$rate.id,
    c(100, 101, 103)
  )

  # explicit id filters
  expect_setequal(
    load_vmdb_rates(con, session.id = c(2))$observations$rate.id,
    c(101, 102)
  )
  expect_equal(
    load_vmdb_rates(con, rate.id = c(102))$observations$rate.id,
    c(102)
  )

  # load_vmdb_magnitudes: without extras
  res <- load_vmdb_magnitudes(con, shower = 'PER')
  expect_type(res, 'list')
  expect_true(is.data.frame(res$observations))
  expect_null(res$sessions)
  expect_equal(nrow(res$observations), 2)

  # with sessions and magnitudes
  res <- load_vmdb_magnitudes(con, shower = 'PER', withSessions = TRUE, withMagnitudes = TRUE)
  expect_true(is.data.frame(res$sessions))
  expect_true(methods::is(res$magnitudes, 'table'))
  expect_equal(sort(row.names(res$sessions)), c('1', '2'))
  expect_equal(sort(row.names(res$magnitudes)), c('200', '201'))

  # magnitude filters mirror rate filters
  expect_equal(
    load_vmdb_magnitudes(con, session.id = c(3))$observations$magn.id,
    c(203)
  )
  expect_equal(
    load_vmdb_magnitudes(con, magn.id = c(202))$observations$magn.id,
    c(202)
  )
  expect_setequal(
    load_vmdb_magnitudes(con, period = c('2015-08-12 00:00:00', '2015-08-12 02:00:00'))$observations$magn.id,
    c(200, 201)
  )
  expect_equal(
    load_vmdb_magnitudes(con, shower = c(NA), sl = c(340.0, 10.0))$observations$magn.id,
    c(203)
  )
})
