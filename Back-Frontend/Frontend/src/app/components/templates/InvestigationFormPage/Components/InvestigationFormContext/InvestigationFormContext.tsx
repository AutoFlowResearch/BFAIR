import { Provider } from 'react-redux';
import store from '../../../../../shared/redux/store';

// const initialState = new Investigation();

// export const InvestigationContext = React.createContext([
//   new Investigation(),
//   (() => {}) as any,
// ]);

const InvestigationFormContextProvider = (props) => {
  return <Provider store={store}>{props.children}</Provider>;
};

export default InvestigationFormContextProvider;
