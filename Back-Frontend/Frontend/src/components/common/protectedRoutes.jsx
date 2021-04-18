import React from "react";
import { Route, Redirect } from "react-router-dom";

export const ProtectedRoute = ({ component: Component, roles, ...rest }) => (
  <Route
    {...rest}
    render={props => {
      // get current user
      const currentUser = localStorage.getItem("token");
      const role = "admin";
      if (!currentUser) {
        // not logged in so redirect to login page with the return url
        return (
          <Redirect
            to={{ pathname: "/login", state: { from: props.location } }}
          />
        );
      }

      // check if route is restricted by role
      if (roles && roles.indexOf(role) === -1) {
        // role not authorised so redirect to home page
        return <Redirect to={{ pathname: "/not-authorized" }} />;
      }

      // authorised so return component
      return <Component {...props} />;
    }}
  />
);
