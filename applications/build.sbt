name := "Cassia"

version := "0.0"

scalaVersion := "2.10.2"

scalaSource in Compile <<= baseDirectory(_ / "src")

scalaSource in Test <<= baseDirectory(_ / "test")

libraryDependencies ++= Seq(
  "org.scala-lang" % "scala-reflect" % "2.10.2"
)

scalacOptions += "-deprecation"

scalacOptions += "-unchecked"

scalacOptions += "-Xlog-free-terms"

fork in run := true

javaOptions in run += "-Djava.library.path=lib/"